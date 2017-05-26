#!/usr/bin/perl -w
use strict;
use Getopt::Long;
# This works as of June 2014
use List::Util qw ( min max sum );
use CXGN::Phylo::Mrbayes;

# read in MrBayes .run*.t files
# remove branch lengths
# put newicks in canonical form (so exactly one newick for each topology)
# store newicks in hash

my $alignment_nex_filename = shift; # e.g. fam4436.nex if trees are in fam4436.nex.run1.t etc.

my $default_burnin_frac = 0.2;
my $start_gen = undef;
my $burnin_frac = undef;
my $n_taxa = undef;
my $n_chars = undef;

GetOptions(
    'burnin=i'           => \$start_gen,
    'burnin_fraction'     => \$burnin_frac, # exclamation point means can use  -nj  and  -nonj
	   'n_taxa=i' => \$n_taxa,
	   'n_chars=i' => \$n_chars, # 'characters' i.e. alignment columns
);

if(!defined $burnin_frac){
  $burnin_frac = $default_burnin_frac;
}

my ($n_runs0, $n_Ts, $delta_T, $n_swaps, $sample_freq0) = (-1, -1, -1, -1, -1); # (2, 4, 0.1, 1, 500);
if(-f $alignment_nex_filename){
open my $fh_nexus, "<", "$alignment_nex_filename";
while(<$fh_nexus>){
$n_runs0 = $1 if(/nruns\s*=\s*(\d+)/);
$n_Ts = $1 if(/nchains\s*=\s*(\d+)/);
$delta_T = $1 if(/temp\s*=\s*([0-9.]+)/);
$n_swaps = $1 if(/nswaps\s*=\s*(\d+)/);
$sample_freq0 = $1 if(/samplefreq\s*=\s*(\d+)/);
};
}else{ # no .nexus file
  die if(!defined $n_chars or !defined $n_taxa);
}



my @runfiles = `ls $alignment_nex_filename.run*.t`;
my $n_runs = scalar @runfiles;
if($n_runs == 0){
print "$alignment_nex_filename  No run*.t files. \n";
exit;
}elsif($n_runs == 1){
  print "$alignment_nex_filename Only one run*.t file. No comparison between runs of topology distribution.\n";
  exit;
}
print "$alignment_nex_filename \n";
warn "Number of runs inconsistency: $n_runs0, $n_runs.\n" if($n_runs != $n_runs0);
my $last_3lines = `tail -3 $runfiles[0]`;
my ($next_to_last_line, $last_line, $really_last_line) = split("\n", $last_3lines);
print "$next_to_last_line \n";
print "$last_line \n";

$last_line =~ /tree gen[.](\d+) =/;
my $last_gen = $1;

$next_to_last_line =~ /tree gen[.](\d+) =/;
my $next_to_last_gen = $1;

my $sample_freq = $last_gen - $next_to_last_gen;
warn "Sample frequency inconsistency: $sample_freq0, $sample_freq. \n" if($sample_freq != $sample_freq0);
if(!defined $start_gen){
  $start_gen = int($burnin_frac * $last_gen);
}

print "# nruns: $n_runs  nchains: $n_Ts  deltaT: $delta_T  nswaps: $n_swaps \n";
print "# run length: $last_gen   burn-in: $start_gen  samplefreq: $sample_freq0 \n";

my $mrb_obj = CXGN::Phylo::Mrbayes->new(
    {
     'alignment_nex_filename' => $alignment_nex_filename,
     'n_runs' => $n_runs, 
     'burnin_fraction' => $burnin_frac,
     'sample_freq' => $sample_freq,
     'n_taxa' => $n_taxa,
     'n_alignment_columns' => $n_chars,
    });

my ($topo_dists, $splits_dists) = $mrb_obj->topo_analyze($start_gen);
printf("# topology. L1: (%6.3f %6.3f %6.3f %6.3f %6.3f) %6.3f  %6s \n" .
       "#         Linf: (%6.3f %6.3f %6.3f %6.3f %6.3f) %6.3f  %6s  \n\n", 
       @$topo_dists[0..5], ($topo_dists->[2] > 1)? $topo_dists->[5] / $topo_dists->[2] : '---',
       @$topo_dists[ 6 .. 11 ], ($topo_dists->[8] > 1)? $topo_dists->[11] / $topo_dists->[8] : '---' );
printf("# splits.   L1: (%6.3f %6.3f %6.3f %6.3f %6.3f) %6.3f  %6s \n" .
       "#         Linf: (%6.3f %6.3f %6.3f %6.3f %6.3f) %6.3f  %6s \n\n", 
       @$splits_dists[ 0 .. 5  ], ($splits_dists->[2] > 1)? $splits_dists->[5] / $splits_dists->[2] : '---',
       @$splits_dists[ 6 .. 11 ], ($splits_dists->[8] > 1)? $splits_dists->[11] / $splits_dists->[8] : '---' );
