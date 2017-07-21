#!/usr/bin/perl
use strict;
use warnings;
#use Getopt::Std;
use Getopt::Long;
use List::Util qw ( min max sum );

#use lib '/home/tomfy/Orthologger/lib';
use lib '/home/tomfy/MB/lib';
use Overlap;
use Mrbayes; # perl module encapsulating mrbayes bayesian phylogeny program.

#use Devel::Cycle; # for finding circular refs - cause of memory leaks.
# find_cycle($test); # to find circular refs in $test (which might be an object, e.g.)

# read in an alignment file. get an overlap object. 

#use vars qw($opt_i $opt_S $opt_f);
# -i <input_alignment_filename>     (name of input alignment file, fasta format)
# -S <seed>  (rng seed, default: get from clock)
# -f <fraction> (fraction of non-gap chars required in alignment column to keep it. Default is 0.8)


# typical usage: perl MB.pl -i fam9877.nex

# get options
#getopt("i:S:f:");

# defaults:
my $fasta_input_file = undef;
my $nex_input_file = undef;
my $seed = srand();
my $swapseed = srand();
my $nongap_fraction = 0.8;
my $chunk_size = 2000;
my $print_freq = 500;
my $n_temperatures = 3;
my $n_temperatures_out = undef; # $n_temperatures;
my $deltaT = 0.1;
my $Tratio = undef;
my $Tladder = undef;
my $sample_freq = 20;
my $n_runs = 2;
my $burnin_fraction = 0.1;
my $converged_chunks_required = 10;
my $modelparam_min_ok_ESS = 200;
my $append = 'yes';
my $max_gens = 1000000;
my $min_chunk_size = 100;       # mb doesn't allow anything less
my $use_mpi = undef;
my $max_processors = 1;
my $mb_name = 'mb';
my $max_ok_L1 = 0.1;
my $n_swaps = 1;    # number of T-swap moves attempted per generation.

my $reproducible = 0; # if true will give identical results each time run, BUT
# the results will be INCORRECT (if multiple chunks). So use reproducible=true ONLY for testing.
GetOptions('fasta_input_file=s' => \$fasta_input_file, # fasta alignment file
	   'seed=i' => \$seed,
           'swapseed=i' => \$swapseed,
	   'nongap_fraction=f' => \$nongap_fraction,
	   'chunk_size=i' => \$chunk_size,
           'print_every=i' => \$print_freq,
	   'nTs=i' => \$n_temperatures,
           'nTs_out=i' => \$n_temperatures_out,
	   'deltaT=s' => \$deltaT,
           'Tratio=f' => \$Tratio,
	   'Tladder=s' => \$Tladder, # e.g. '1.0, 1.2, 1.5';
	   'sample_freq=i' => \$sample_freq,
	   'n_runs=i' => \$n_runs,
	   'burn-in_fraction=f' => \$burnin_fraction,
	   'converged_chunks_required=i' => \$converged_chunks_required,
	   'ESS_min=i' => \$modelparam_min_ok_ESS,
           'append=s' => \$append,
           'max_gens=i' => \$max_gens,
           'reproducible=i' => \$reproducible,
	   'use_mpi=i' => \$use_mpi,
           'max_processors=i' => \$max_processors,
           'mb_name=s' => \$mb_name,
           'n_swap=i' => \$n_swaps,
           'max_ok_L1=f' => \$max_ok_L1,
           'nex_input_file=s' => \$nex_input_file,
          );

if ($chunk_size < $min_chunk_size) {
   warn "Resetting chunk size form $chunk_size to $min_chunk_size (min allowed by MrBayes)\n";
   $chunk_size = $min_chunk_size;
}


if (defined $Tladder) {
   $Tladder =~ s/\s//g;         # remove whitespace
   # otherwise leave as is 
   my @Ts = split(",", $Tladder);
   $n_temperatures = scalar @Ts; # if Tladder specified, get n_temperatures from it.
} elsif (defined $Tratio) {
   my $T = 1.0;
   my @Ts = ();
   for (1..$n_temperatures) {
      push @Ts, $T;
      $T *= $Tratio;
   }
   $Tladder = '(' . join(',', @Ts) . ')';
} else {
   my $T = 1.0;
   my @Ts = ();
   for (1..$n_temperatures) {
      push @Ts, $T;
      $T += $deltaT;
   }
   $Tladder = '(' . join(',', @Ts) . ')';
}

$n_temperatures_out = $n_temperatures if(!defined $n_temperatures_out); # default is output all temperatures

print "Seed: $seed\n";
print "chunksize: $chunk_size \n";
print "nongapfrac: $nongap_fraction\n";
print "n_temperatures: $n_temperatures\n";
print "Tladder: $Tladder \n";
print "sample_freq: $sample_freq\n";
print "n_runs: $n_runs\n";
print "burnin fraction: $burnin_fraction\n";
print "converged chunks required: $converged_chunks_required\n";
print "modelparam min OK ESS: $modelparam_min_ok_ESS\n";
print "append: $append \n";
print "max processors: $max_processors \n";
# exit;

die "Must specify name of input alignment file. Either -nex_input [nexus_filename], or -fasta_input [fast_filename].\n" . "Usage: MB.pl --nex_input_file <nex_input_filename> [--seed <rng seed> --nongap_fraction <overlap fraction>]\n" unless(($nex_input_file and -f $nex_input_file) or ($fasta_input_file and -f $fasta_input_file));

#my $input_file = $opt_i;




# construct MrBayes object and run
#my $seed = ($opt_S)? $opt_S : undef;
#my $swapseed = ($seed)? $seed + 1000 : undef;
#my $n_temperatures = 4;
# $n_swaps = ($n_temperatures > 2)? int ( ($n_temperatures-1) * $n_temperatures/2 / 2.1 ) : 1; # 1,2,3 T's -> 1, 4 T's ->2, 5 T's -> 4
print "n_swaps (number of swap attempts per generation): $n_swaps \n";
#my $n_runs = 2;
#my $temperature_gap = 0.3;   # temperature spacing of chains
#my $chunk_size = 200;	   # steps between decisions whether to go on.
#my $print_freq = int($chunk_size/2);
#my $sample_freq = 20; # write out state every this many steps.
$sample_freq = max($sample_freq, 1);


my $alignment_nex_filename;
if (defined $nex_input_file) {
   $alignment_nex_filename = $nex_input_file;
} else {
   my $mrb_outfile_basename;
   $mrb_outfile_basename = $fasta_input_file;
   $mrb_outfile_basename =~ s/[.]fasta//;
   $alignment_nex_filename = $mrb_outfile_basename . '.nex';
   #### Get the alignment:
   my $align_string =  `cat $fasta_input_file`;
   # fixes to $align_string:
   $align_string =~ s/IMGA[|]/IMGA_/g; #pipes in id cause problem; replace '|' with '_'.
   my $fixprefix = 'X_'; # if id begins with non-alphabetic char, prefix with this.
   $align_string =~ s/^>(\s*)([^a-zA-Z])/>$1$fixprefix$2/xmsg; # to make clearcut happy.

   # construct an overlap object.
   # my $bootstrap_seed = 1234567;	# ($opt_S)? $opt_S : undef;
   #my $nongap_fraction = ($opt_f)? $opt_f : 0.8;
   # print "$align_string \n";
   my $overlap_obj = Overlap->new($align_string, $nongap_fraction); # , $bootstrap_seed);
 
   my $overlap_nexus_string = $overlap_obj->overlap_nexus_string();

   # print $overlap_nexus_string, "\n";
   open my $fh1, ">", "$alignment_nex_filename";
   print $fh1 $overlap_nexus_string, "\n";
   close $fh1;
}


print "Seed, swapseed: $seed  $swapseed \n"; #sleep(1);
# print "reproducible: $reproducible \n";
my $mrb_obj = Mrbayes->new(
			   {'alignment_nex_filename' =>$alignment_nex_filename,
			    'chunk_size' => $chunk_size,
                            'print_freq' => $print_freq,
			    'seed' => $seed,
			    'swapseed' => $swapseed,
			    #	   'fixed_pinvar' => 0.4
			    'sample_freq' => $sample_freq,
			    'modelparam_min_ok_ESS' => $modelparam_min_ok_ESS,
			    'n_temperatures' => $n_temperatures,
                            'n_temperatures_out' => $n_temperatures_out,
			    'temperatures_parameter' => $Tladder,
                            #		    'user_defined_temperatures' => $user_def_temperatures,
			    'n_swaps' => $n_swaps,
			    'n_runs' => $n_runs,
			    'burnin_fraction' => $burnin_fraction,
			    'converged_chunks_required' => $converged_chunks_required,

			    'append' => $append,
			    'max_gens' => $max_gens,
			    'reproducible' => $reproducible,
			    'use_mpi' => $use_mpi,
			    'max_processors' => $max_processors,
			    'mb_name' => $mb_name,
			    'modelparam_min_ok_ESS'     => 2000,
                            'modelparam_max_ok_KSD' => 0.2,
                            'max_ok_L1'             => $max_ok_L1, #    0.01,

                            #			    'temperature_factor' => 1.414, # I wanted to have T's exponentially spaced, but mb does not allow
			   }
			  );
print STDERR "Now mrb_obj->run() \n";
$mrb_obj->run();

exit;
