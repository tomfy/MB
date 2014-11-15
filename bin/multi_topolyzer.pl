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
#
#my $alignment_nex_filename = shift;           # e.g. fam4436.nex if trees are in fam4436.nex.run1.t etc.

my $chunk_topology_count = {};	# hashref of chunk-number;hashref 
my $topology_count = {};	# hashref of chunk-number;hashref 

my $default_burnin_frac = 0.2;
my $start_gen = undef;
my $burnin_frac = undef;
my $prefix = 'chunk';
my ($begin_chunk_number, $end_chunk_number) = (1, 1000);
GetOptions(
	   'burnin=i'           => \$start_gen,
	   'burnin_fraction=f'     => \$burnin_frac, # exclamation point means can use  -nj  and  -nonj
	   'prefix=s' => \$prefix,
	   'begin=i' => \$begin_chunk_number,
	   'end=i' => \$end_chunk_number,
	  );

if (!defined $burnin_frac) {
  $burnin_frac = $default_burnin_frac;
}



for my $chunk_number ($begin_chunk_number .. $end_chunk_number) {
  my $alignment_nex_filename = "$prefix$chunk_number.nexus";
  next if(! -f $alignment_nex_filename);

  my ($n_runs0, $n_Ts, $delta_T, $n_swaps, $sample_freq0) = (2, 4, 0.1, 1, 500);
  open my $fh_nexus, "<", "$alignment_nex_filename";
  while (<$fh_nexus>) {
    $n_runs0 = $1 if(/nruns\s*=\s*(\d+)/);
    $n_Ts = $1 if(/nchains\s*=\s*(\d+)/);
    $delta_T = $1 if(/temp\s*=\s*([0-9.]+)/);
    $n_swaps = $1 if(/nswaps\s*=\s*(\d+)/);
    $sample_freq0 = $1 if(/samplefreq\s*=\s*(\d+)/);
  }
  ;


  my @runfiles = `ls $alignment_nex_filename.run*.p`;
  my $n_runs = scalar @runfiles;
  warn "Number of runs inconsistency: $n_runs0, $n_runs.\n" if($n_runs != $n_runs0);
  my $last_2lines = `tail -2 $runfiles[0]`;
  my ($next_to_last_line, $last_line) = split("\n", $last_2lines);
  my @cols = split(" ", $last_line);
  my $last_gen = $cols[0];
  @cols = split(" ", $next_to_last_line);
  my $next_to_last_gen = $cols[0];
  my $sample_freq = $last_gen - $next_to_last_gen;
  warn "Sample frequency inconsistency: $sample_freq0, $sample_freq. \n" if($sample_freq != $sample_freq0);
  if (!defined $start_gen) {
    $start_gen = int($burnin_frac * $last_gen);
  }


  print "# nruns: $n_runs0  nchains: $n_Ts  deltaT: $delta_T  nswaps: $n_swaps \n";
  print "# run length: $last_gen   burn-in: $start_gen  samplefreq: $sample_freq0 \n";

  my $mrb_obj = CXGN::Phylo::Mrbayes->new(
					  {
					   'alignment_nex_filename' => $alignment_nex_filename,
					   'n_runs' => $n_runs, 
					   'burnin_fraction' => $burnin_frac,
					   'sample_freq' => $sample_freq,
					  });

  my ($topo_dists, $splits_dists) = $mrb_obj->topo_analyze($start_gen);
  # printf("# topology. L1: (%6.3f %6.3f %6.3f %6.3f %6.3f) %6.3f  %6s \n" .
  # 	 "#         Linf: (%6.3f %6.3f %6.3f %6.3f %6.3f) %6.3f  %6s  \n\n", 
  # 	 @$topo_dists[0..5], ($topo_dists->[2] > 1)? $topo_dists->[5] / $topo_dists->[2] : '---',
  # 	 @$topo_dists[ 6 .. 11 ], ($topo_dists->[8] > 1)? $topo_dists->[11] / $topo_dists->[8] : '---' );
  # printf("# splits.   L1: (%6.3f %6.3f %6.3f %6.3f %6.3f) %6.3f  %6s \n" .
  # 	 "#         Linf: (%6.3f %6.3f %6.3f %6.3f %6.3f) %6.3f  %6s \n\n", 
  # 	 @$splits_dists[ 0 .. 5  ], ($splits_dists->[2] > 1)? $splits_dists->[5] / $splits_dists->[2] : '---',
  # 	 @$splits_dists[ 6 .. 11 ], ($splits_dists->[8] > 1)? $splits_dists->[11] / $splits_dists->[8] : '---' );
  my $this_chunk_topology_count = $mrb_obj->get_topology_count();
  $chunk_topology_count->{$chunk_number} = $this_chunk_topology_count;
  while (my ($topo, $count) = each %$this_chunk_topology_count) {
    $topology_count->{$topo} += $count;
  }

}
my $filename = "topology_pd_vs_chunks_$begin_chunk_number-$end_chunk_number";
open my $fhout, ">", "$filename";
print $fhout "# Topology posterior distributions for chunks $begin_chunk_number through $end_chunk_number.\n", 
  "# Distinct topologies found: ", scalar keys %$topology_count, "\n";
my @skeys = sort 
  # { $topology_count->{$b} <=> $topology_count->{$a} } 
  { $b cmp $a}
  keys %$topology_count;
for (@skeys) {
  print $fhout "# $_   ", $topology_count->{$_}, "\n";
}

my @chunk_numbers = sort {$a <=> $b } keys %$chunk_topology_count;

for my $chunk_number (@chunk_numbers) {
  my $topo_count = $chunk_topology_count->{$chunk_number};

  printf $fhout ("%4i    ", $chunk_number);
my @topo_counts = map( (exists $topo_count->{$_})? $topo_count->{$_} : 0, @skeys );

  my $sum_count_over_topos = sum(@topo_counts);

  for (@skeys) {
    my $count = (exists $topo_count->{$_})? $topo_count->{$_} : 0;
    printf $fhout ("%6i ", $count);
#$sum_count_over_topos += $count;
  }
  print $fhout "  $sum_count_over_topos \n";
}
close $fhout;



