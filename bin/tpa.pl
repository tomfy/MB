#!/usr/bin/perl
use strict;
use warnings;
#use Getopt::Std;
use Getopt::Long;
use List::Util qw ( min max sum );
use File::Spec qw ( splitpath );
#use lib '/home/tomfy/Orthologger/lib';
use lib '/home/tomfy/MB/lib';
#use Overlap;
use Mrbayes; # perl module encapsulating mrbayes bayesian phylogeny program.

my %orderednewick_count = ();
my %newick_count = ();

my %newick_number_map  = ();    # %{ $self->{newick_number_map} };
my %number_newick_map  = ();    # %{ $self->{number_newick_map} };

my %number_splits_map = ();

open my $fhxxx, ">", "XXXXX"; 

my $tp_filename = undef;
my $data_nex_filename = undef;
my $burn_in = 0.1;
my $output_factor = 1.4;
GetOptions(
           'in_file=s' => \$tp_filename, # a .tp file or .otps file 
           'nexus_filename=s' => \$data_nex_filename,
           'burn_in=f' => \$burn_in,
           'output_factor=f' => \$output_factor,
          );

#my $leaf_count = n_leaves($data_nex_filename);

my ($n_runs, $n_temperatures, $sample_spacing, $out_lines_per_out_generation) = get_nruns_nTs($tp_filename);
print "# $n_runs  $n_temperatures   $sample_spacing  $out_lines_per_out_generation\n";

my $n_temperatures_out = $n_temperatures;

my ($all_lines, $last_gen) = read_tp_file($tp_filename);
# get number of leaves from newick expression;
my $first_line = $all_lines->[0];
my @cols = split(" ", $first_line);
my $nwck = $cols[9];
if($nwck =~ /.+:(.*):.+/){
   $nwck = $1;
}
my $leaf_count = 0;
while ($nwck =~ s/\d+//) {
   $leaf_count++;
}
print "# leaf count: $leaf_count\n";

my $max_burn_in_gen = ($burn_in < 1.0)? int($burn_in*$last_gen) : $burn_in;
$max_burn_in_gen = $max_burn_in_gen - ($max_burn_in_gen % $sample_spacing);
print "# burn-in: $max_burn_in_gen \n";

my $Topology_chain_data_u = ChainData->new( # unweighted
                                           {
                                            'parameter_name' => 'Topology',
                                            'n_runs'         => $n_runs,
                                            'n_temperatures' => $n_temperatures,
                                            'n_temperatures_out' => $n_temperatures_out,
                                            'gen_spacing'    => $sample_spacing,
                                            'binnable'       => 0
                                           }
                                          ); #
my $Topology_chain_data_w = ChainData->new(  # weighted
                                           {
                                            'parameter_name' => 'Topology',
                                            'n_runs'         => $n_runs,
                                            'n_temperatures' => $n_temperatures,
                                            'n_temperatures_out' => $n_temperatures_out,
                                            'gen_spacing'    => $sample_spacing,
                                            'binnable'       => 0
                                           }
                                          ); #


my $Splits_chain_data_u = ChainData->new( # unweighted
                                         {
                                          'parameter_name' => 'Splits',
                                          'n_runs'         => $n_runs,
                                          'n_temperatures' => $n_temperatures,
                                          'n_temperatures_out' => $n_temperatures_out,
                                          'gen_spacing'    => $sample_spacing,
                                          'binnable'       => 0,
                                          'set_size'       => $leaf_count - 3 # n interior edges, i.e. number of splits for each topology
                                         }
                                        );
my $Splits_chain_data_w = ChainData->new( # weighted
                                         {
                                          'parameter_name' => 'Splits',
                                          'n_runs'         => $n_runs,
                                          'n_temperatures' => $n_temperatures,
                                          'n_temperatures_out' => $n_temperatures_out,
                                          'gen_spacing'    => $sample_spacing,
                                          'binnable'       => 0,
                                          'set_size'       => $leaf_count - 3 # n interior edges, i.e. number of splits for each topology
                                         }
                                        );

# get some ISWs from close to the end of run:
my $maxlnISW =  get_lnISWs($all_lines, $n_temperatures);
for my $it (1 .. scalar keys %$maxlnISW) {
   print "$it  ", $maxlnISW->{$it-1}, "  ";
}
print "\n";

my $topology_number = 0;
my $n_topologies_so_far = 0;
my $Topology_count = {};
my $end_gen = 2*$max_burn_in_gen;

for my $aline (@$all_lines) {
   my @cols = split(" ", $aline);
   my ($gen, $i_run, $i_T, $i_w, $lnISW, $newick, $the_splits);
   if (scalar @cols == 11) {    # ordered, long format
      ($gen, $i_run, $i_T, $i_w, $lnISW, $newick, $the_splits) = @cols[0,5,6,7,8,9,10];
      if ( !exists $newick_number_map{$newick} ) { # if a new topology, get its number and store newick<->number maps
         $n_topologies_so_far++; # counts distinct topologies
         $newick_number_map{$newick}          = $n_topologies_so_far; # 1,2,3,...
         $number_newick_map{$n_topologies_so_far} = $newick;
         #  print "Tnumber, newick: $topology_number $newick \n";
         $topology_number = $n_topologies_so_far;
      } else {
         $topology_number = $newick_number_map{$newick};
      }
   } elsif (scalar @cols == 10) { # raw, or ordered, short format
      if ($cols[9] =~ /[:]/) {    # ordered, short format 
         ($gen, $i_run, $i_T, $i_w, $lnISW) =  @cols[0,5,6,7,8];
         ($topology_number, $newick, $the_splits) = split(":", $cols[9]);
         if ($newick ne '-') {
            $number_newick_map{$topology_number} = $newick;
       #     $newick_number_map{$newick} = $topology_number;
            $number_splits_map{$topology_number} = $the_splits;
         } else {
            $newick = $number_newick_map{$topology_number};
            $the_splits = $number_splits_map{$topology_number};
         }
      #   print "ABC:  $topology_number $newick $the_splits \n";
      } else {                  # newick is raw, not yet ordered
         ($gen, $i_run, $i_T, $i_w, $lnISW, $newick) = @cols[0,5,6,7,8,9];
         #    print "$gen $newick $leaf_count\n";
         my ($minlabel, $ordered_newick, $tree_labels, $split_count) = get_ordered_newick($newick, $leaf_count);
         $newick = $ordered_newick;
         $newick =~ s/[']//g;
         $the_splits = join(";", sort {$a <=> $b} keys %$split_count);
         #    print "xxx $the_splits \n";
         if ( !exists $newick_number_map{$newick} ) { # if a new topology, get its number and store newick<->number maps
            $n_topologies_so_far++; # counts distinct topologies
            $newick_number_map{$newick}          = $n_topologies_so_far; # 1,2,3,...
            $number_newick_map{$n_topologies_so_far} = $newick;
            #  print "Tnumber, newick: $topology_number $newick \n";
            $topology_number = $n_topologies_so_far;
         } else {
            $topology_number = $newick_number_map{$newick};
         }
      }
   } else {
      die "unknown format.\n";
   }
   #print $fhxxx join(' ', @cols[0..8]), "    ", $newick, "    ", $the_splits, "\n";

   $i_T %= $n_temperatures;

   next if($cols[0] <= $max_burn_in_gen); # skip the burn-in

   # store data points:
   my $ISW = 1000*exp($lnISW - $maxlnISW->{$i_T});
   my $sample_set_id = $i_run . "_" . ($i_T) . "_" . ($i_w % $n_temperatures); # walkers separate

   my $ntopo = $topology_number; #$newick_number_map{$newick};
   $Topology_chain_data_u
     ->store_data_point( $sample_set_id, $gen, $ntopo);
 #  print "ntopo: $ntopo\n";
 #  print "newick number: ", $newick_number_map{$newick}, "\n";
   my $ntopo_w = $ntopo . ":$ISW"; # this is 
#   print "$gen  $i_run $i_T $i_w   Ntopo, ISW: $ntopo  ", $newick_number_map{$newick}, "  $ISW  \n" if($ntopo != $newick_number_map{$newick});
   $Topology_chain_data_w
     ->store_data_point( $sample_set_id, $gen, $ntopo_w); # $newick_number_map{$newick} );
   $Topology_count->{$newick}++;

   # store split
   $Splits_chain_data_u->store_data_point($sample_set_id, $gen, $the_splits );
   $Splits_chain_data_w->store_data_point($sample_set_id, $gen, $the_splits . ":$ISW" );

   #  printf $fhotps ("%6d  %s %s %s %s  %3d %3d %3d   %10.3f  %s  %s\n", @cols[0..8], $ordered_newick, $the_split);
   #print "$gen $end_gen \n";
   # ****************** ANALYZE ********************
   if (($gen == $end_gen)  and  ($i_run == $n_runs-1) and ($i_T == $n_temperatures_out-1)) { # the last line of this generation
      # analyze data up to this point and output

      open my $fhhist, ">", "Topology_histogram_u";
      print $fhhist $Topology_chain_data_u->{histograms}->histogram_string('by_bin_weight', 1.0) if($gen == $last_gen);
      #    my @topo_L1_distances_u = $Topology_chain_data_u->{histograms}->avg_L1_distance();
      #    printf $fhhist ( "# L1:   %6.3f %6.3f %6.3f %6.3f %6.3f\n# Linf: %6.3f %6.3f %6.3f %6.3f %6.3f \n\n", @topo_L1_distances_u[ 0 .. 4 ], @topo_L1_distances_u[6..10] );
      close $fhhist;

      open $fhhist, ">", "Topology_histogram_w";
      print $fhhist $Topology_chain_data_w->{histograms}->histogram_string('by_bin_weight', 1.0) if($gen == $last_gen);
      #   my @topo_L1_distances_w = $Topology_chain_data_w->{histograms}->avg_L1_distance();
      #   printf $fhhist ( "# L1:   %6.3f %6.3f %6.3f %6.3f %6.3f\n# Linf: %6.3f %6.3f %6.3f %6.3f %6.3f \n\n", @topo_L1_distances_w[ 0 .. 4 ], @topo_L1_distances_w[6..10] );
      close $fhhist;

      #   printf ( " gen: %6d   L1u:   %6.3f %6.3f %6.3f %6.3f %6.3f    L1w: %6.3f %6.3f %6.3f %6.3f %6.3f \n", $gen, @topo_L1_distances_u[ 0 .. 4 ], @topo_L1_distances_w[0..4] );

      my ($l1s, $l2s) = $Topology_chain_data_w->{histograms}->get_l12s();
      my ($mean_l1, $var_l1, $stderr_l1) = mean_variance_stderr($l1s);
      my ($mean_l2, $var_l2, $stderr_l2) = mean_variance_stderr($l2s);
      #print "l1s: ", join("  ", @$l1s), "  l2s: ", join("  ", @$l2s), "\n";
   #   printf join(" ", @$l1s), "  ";
      printf("gen: %6d   l1: %6.4f +- %6.4f;  l2: %7.5f +- %7.5f   ",
             #   %6d %6d %6d %6d \n", 
             $gen, $mean_l1, $stderr_l1, $mean_l2, $stderr_l2);
      # scalar keys %newick_count, sum(values %newick_count), scalar keys %orderednewick_count, sum(values %orderednewick_count));

      open $fhhist, ">", "Splits_histogram_u";
      print $fhhist $Splits_chain_data_u->{histograms}->histogram_string('by_bin_weight', 1.0);
      #   my @splits_L1_distances_u = $Splits_chain_data_u->{histograms}->avg_L1_distance();
      #   printf $fhhist ("# L1:   min: %5.3f; q1: %5.3f; median: %5.3f; q3: %6.3f; max %5.3f \n",  @splits_L1_distances_u[ 0 .. 4 ]);
      #   printf $fhhist ("# Linf: min: %5.3f; q1: %5.3f; median: %5.3f; q3: %5.3f; max: %5.3f \n\n", @splits_L1_distances_u[ 6 .. 10 ] );
      close $fhhist;

      open $fhhist, ">", "Splits_histogram_w";
      print $fhhist $Splits_chain_data_w->{histograms}->histogram_string('by_bin_weight', 1.0);
      #  my @splits_L1_distances_w = $Splits_chain_data_w->{histograms}->avg_L1_distance();
      #   printf $fhhist ("# L1:   min: %5.3f; q1: %5.3f; median: %5.3f; q3: %6.3f; max %5.3f \n",  @splits_L1_distances_w[ 0 .. 4 ]);
      #   printf $fhhist ("# Linf: min: %5.3f; q1: %5.3f; median: %5.3f; q3: %5.3f; max: %5.3f \n\n", @splits_L1_distances_w[ 6 .. 10 ] );
      close $fhhist;

      my ($rravgd1s, $w_avg_d1s) = compare_runs_each_T($Topology_chain_data_w, $n_runs, $n_temperatures);
      while (my ($it, $d1) = each @$rravgd1s) {
         printf("%3d %7.5f  ", $it, $d1);
      }
      print "  ";
      while (my ($it, $d1) = each @$w_avg_d1s) {
         printf("%3d %7.5f  ", $it, $d1);
      }

      print "  ", compare_allT_estimates($Topology_chain_data_w, $w_avg_d1s), "\n";

      # get next analysis gen:
      $end_gen = min(int($output_factor*$end_gen) + $sample_spacing, $last_gen);
      $end_gen -= $end_gen % $sample_spacing;
   }
}
open my $fhtrfj, ">", 'TRFjumps';
print $fhtrfj analyze_chain_splits($Splits_chain_data_u->{setid_gen_value}, $n_runs, $n_temperatures);

sub TRF_distance{ # Topological Robinson-Foulds distance
   # $spl1, $spl2 are  ':' separated lists of bit patterns indicating presence/absence of leaves in split (always the odd one)
   my $spl1 = shift;            # string, e.g.:  5:253:1393:2921 
   my $spl2 = shift;
   my @splits1 = split(";", $spl1);
   my @splits2 = split(";", $spl2);
   my %s1 = map(($_ => 1), @splits1);
   my %s2 = map(($_ => 1), @splits2);
# print STDERR "XXX:  ", join(";", keys %s1), "   ", join(";", keys %s2), "\n";

   my $one_only_count = 0;
   my $both_count = 0;
   for my $the_split (keys %s1) {
      if (exists $s2{$the_split}) {
         $both_count++;
      } else {
         $one_only_count++;
      }
   }
#   print "both_count: $both_count  only_one_count: $one_only_count \n";
   return $one_only_count;
}

# my ($topology_chain_data, $splits_chain_data, $number_newick_map) = 
#  retrieve_topology_samples_from_tp($tp_filename, $leaf_count, $n_runs, $n_temperatures, $n_temperatures_out, $sample_spacing, $burn_in_fraction);
#print "tp_filename: $tp_filename\n";
# tp_analyze($topology_chain_data, $splits_chain_data, $number_newick_map, $tp_filename, $n_runs, $n_temperatures, $n_temperatures_out);
#$mbobj->retrieve_topology_samples_from_tp( $tp_filename, 0);


# print "distinct topologies ", scalar keys %orderednewick_count, "  ", scalar keys %newick_count, "\n";
# my $dot_tp_filename = shift; # e.g. fam9877.nex
#  print STDERR "dot_tp_filename: $dot_tp_filename \n";
#  my $n_taxa = shift;
#  my $n_runs = shift;             #= $self->{n_runs};
#  my $n_temperatures = shift;     #$self->{n_temperatures};
#  my $n_temperatures_out = shift; #$self->{n_temperatures_out};
#  my $sample_freq = shift;
#  my $burn_in_fraction = shift // 0.1;

sub compare_runs_each_T{ 
# 1) lump together walkers within a temperature, then for each temperature, compare each run to mean of all runs.
# 2) compare all n_T * n_runs walkers at each temperatures, compare each walker to mean of all.
   my $the_chain_data = shift;
   my $n_runs = shift;
   my $n_temperatures = shift;
   my @d1s1 = ();
   my @d1s2 = ();
   for my $Tindex (0..$n_temperatures_out-1) { # lump together the walkers of each run & temperature.
         my ($compare_string1, $compare_string2) = ('', '');
         my @grp_strs1 = ();
      #   my @grp_strs2 = ();
         my @rtws2 = ();
         for my $Rindex (0..$n_runs-1) {
            my $rt_bit = $Rindex . '_' . $Tindex . '_';
            my @rtws1 = ();
            for my $Windex (0..$n_temperatures-1) {
               my $rtw_bit = $rt_bit . $Windex;
               push @rtws1, $rtw_bit;
               push @rtws2, $rtw_bit;
            }
            push @grp_strs1, join(',', @rtws1);
         }
         $compare_string1 = join(';', @grp_strs1);
         #print "compare string: $compare_string\n";
         my ($d1_avg, $d2_avg, $dinf_avg) = $the_chain_data->{histograms}->Ln_distances_y($compare_string1);
         push @d1s1, $d1_avg;

         $compare_string2 = join(';', @rtws2);
      #   print "compare string 2: $compare_string2 \n";
  ($d1_avg, $d2_avg, $dinf_avg) = $the_chain_data->{histograms}->Ln_distances_y($compare_string2);
         push @d1s2, $d1_avg;
      }
   return (\@d1s1, \@d1s2);
}

sub compare_allT_estimates{
   my $the_chain_data_w = shift;
   my $avg_w_v_mean_d1s = shift;
   my $n_runs = $the_chain_data_w->{n_runs};
   my $n_Ts =  $the_chain_data_w->{n_temperatures};

   my $weights_string = '';
   my @grp_strs = ();
   for my $Rindex (0..$n_runs-1) {
      my @rtws = ();
      for my $Tindex (0..$n_Ts-1) {
         my $T_weight = 1/$avg_w_v_mean_d1s->[$Tindex]**2;
         $weights_string .= sprintf("%5d   %8.4g  ", $Tindex,  $T_weight) if($Rindex == 0);
         for my $Windex (0..$n_Ts-1) {
            my $rtw_bit = $Rindex . '_' . $Tindex . '_' . $Windex . ':' . $T_weight;
            push @rtws, $rtw_bit;
         }
      }
      push @grp_strs, join(",", @rtws);
   }
#   print "$weights_string \n";
   my $compare_string = join(";", @grp_strs);
   my ($d1_avg, $d2_avg, $dinf_avg) = $the_chain_data_w->{histograms}->Ln_distances_y($compare_string);
   return $d1_avg;              # , $d2_avg, $dinf_avg)
}

sub get_lnISWs{
   my $all_lines = shift;
   my $n_temperatures = shift;
   my %maxlnISW = ();
   #for (@$all_lines[-$out_lines_per_out_generation .. -1]) {
   my $i = -1;
   my $n_lines = scalar @$all_lines;
   while ($i > 1 - $n_lines) {
#  while(1){
      my $a_line = $all_lines->[$i];
      my @cols = split(" ", $a_line);
      my ($i_T, $lnISW) = @cols[6,8];
      $i_T %= $n_temperatures;
      $maxlnISW{$i_T} = (exists $maxlnISW{$i_T})? max($lnISW, $maxlnISW{$i_T}) : $lnISW;
 #    last if($i_T == 0);
      $i--;
   }
   return \%maxlnISW;
}

sub get_ordered_newick{
   my $newick = shift;
   my $n_taxa = shift;
   $newick =~ s/;\s*$//;        # remove final semi-colon.
   $newick =~ s/^\s+//;         # remove init whitespace.
   $newick =~ s/:[0-9]+[.][0-9]+(e[-+][0-9]{2,3})?(,|[)])/$2/g; # remove branch lengths
   my $mask = ( 1 << $n_taxa ) - 1; # n_taxa 1's
   my $split_count = {}; # keys are split bit patterns, values: counts
   my ( $minlabel, $ordered_newick, $tree_labels, $split_bp ) =
     Mrbayes::order_newick( $newick, $split_count, $mask, 0 );
   $newick = $ordered_newick;
   if ( $ordered_newick =~ s/^\s*[(]1,/'(1,(/ ) { # go from form (1,(...),(...)) with trifurcation, to (1,((...),(...))) with only bifurcations
      $ordered_newick .= ")'";
   } else {
      die "newick: $ordered_newick doesn't start with '(1,' as expected.\n";
   }
   return ($minlabel, $ordered_newick, $tree_labels, $split_count);
}

sub get_nruns_nTs{
   my $tp_filename = shift;
   open my $fh, "<", "$tp_filename" or die "Couldn't open $tp_filename for reading.\n";
   my %rtw_count = ();
   my %r_count = ();
   my %t_count = ();
   my %w_count = ();
   my $out_lines_per_generation_out = 0; # number of lines per generation output (should be $n_runs*$n_T_out)
   my $n_lines_read = 0;
   my $old_gen = 0;
   my $sample_spacing = 0;
   my $prev_gen = undef;
   while (<$fh>) {
      next if(/^\s*#/);
      $n_lines_read++;
      my @cols = split(" ", $_);
      my ($gen, $i_run, $i_T, $i_w) = @cols[0, 5,6,7];
      if (defined $prev_gen and $gen != $prev_gen) {
         $sample_spacing = max($gen - $prev_gen, $sample_spacing);
         last;
      }
      $prev_gen = $gen;
      $out_lines_per_out_generation++;
      #print "sss $sample_spacing \n";
      $r_count{$i_run}++;
      $t_count{$i_T}++;
      $w_count{$i_w}++;
      my $rtw = $i_run . '_' . $i_T . '_' . $i_w;
      # $rtw_count{$rtw}++;
      # last if($n_lines_read >= 100  and  $rtw_count{$rtw} == 2);
   }
   close $fh;

   my $n_runs = scalar keys %r_count;
   my $max_T_index = -1;
   for (keys %t_count) {
      $max_T_index = max($max_T_index, int($_ / $n_runs));
   }
   # my $max_w_index = -1;
   # for (keys %w_count) {
   #    $max_w_index = max($max_w_index, int($_ / $n_runs));
   # }

   return ($n_runs, $max_T_index+1, $sample_spacing, $out_lines_per_out_generation);
}

sub n_leaves{
   my $filename = shift;
   open my $fh, "<", "$filename" or die "couldn't open $filename for reading.\n";
   while (<$fh>) {
      if (/dimensions\s+ntax=(\d+)/) {
         return $1;
      }
   }
   return undef;
}

sub read_tp_file{
   my $tp_filename = shift;
   #   my $burn_in = shift || 0.1; # if < 1 interpret as burn-in fraction, else as max burn-in generation number.
   open my $fht, "<", "$tp_filename" or die "Couldn't open $tp_filename for reading.\n";
   my @all_lines = <$fht>;
   my $last_line = $all_lines[-1];
   my $last_gen =  ($last_line =~ /^\s*(\d+)/)? $1 : undef;
   #   my $max_burn_in_gen = ($burn_in < 1)? int($burn_in * $1) : $burn_in;
   return (\@all_lines, $last_gen);
}


sub mean_variance_stderr{
   my $x = shift;               # array ref of numbers.
   my ($count, $sum_x, $sum_xsqr) = (scalar @$x, 0, 0);
   return (undef, undef, undef) if($count == 0);
   for (@$x) {
      $sum_x += $_;
      $sum_xsqr += $_*$_;
   }
   my $mean = $sum_x/$count;
   my $variance = ($count > 1)? ($count/($count-1)) * $sum_xsqr/$count - $mean**2 : undef; # 'bessel correction n/(n-1) applied -> unbiased estimator of population variance.
   my $stderr = ($count > 1)? sqrt($variance/$count) : undef;
   return ($mean, $variance, $stderr);
}


# sub retrieve_topology_samples_from_tp {

#    # read data from ... .tp file
#    # store each param data in separate array of gen/paramval hashrefs
#    my $dot_tp_filename = shift; # e.g. fam9877.nex
#    print STDERR "dot_tp_filename: $dot_tp_filename \n";
#    my $n_taxa = shift;
#    my $n_runs = shift;             #= $self->{n_runs};
#    my $n_temperatures = shift;     #$self->{n_temperatures};
#    my $n_temperatures_out = shift; #$self->{n_temperatures_out};
#    my $sample_freq = shift;
#    my $burn_in_fraction = shift // 0.1;
#    my %newick_number_map  = (); # %{ $self->{newick_number_map} };
#    my %number_newick_map  = (); # %{ $self->{number_newick_map} };
#    my $topology_number    = 0;
#    my $generation;
#    my $newick;
#    my @generations = ();
#    #my %newick_count = ();
#    my $Topology_count = {};

#    # Read param values from .t file and store in ChainData objects.
#    my $Topology_chain_data = ChainData->new(
#                                             {
#                                              'parameter_name' => 'Topology',
#                                              'n_runs'         => $n_runs,
#                                              'n_temperatures' => $n_temperatures,
#                                              'n_temperatures_out' => $n_temperatures_out,
#                                              'gen_spacing'    => $sample_freq,
#                                              'burn_in_fraction' => $burn_in_fraction,
#                                              'binnable'       => 0
#                                             }
#                                            ); #

#    my $Splits_chain_data = ChainData->new(
#                                           {
#                                            'parameter_name' => 'Splits',
#                                            'n_runs'         => $n_runs,
#                                            'n_temperatures' => $n_temperatures,
#                                            'n_temperatures_out' => $n_temperatures_out,
#                                            'gen_spacing'    => $sample_freq,
#                                            'burn_in_fraction' => $burn_in_fraction,
#                                            'binnable'       => 0,
#                                            'set_size'       => $n_taxa - 3 # n interior edges, i.e. number of splits for each topology
#                                           }
#                                          );




#    open my $fht, "<", "$dot_tp_filename" or die "Couldn't open $dot_tp_filename for reading.\n";
#    my @all_lines = <$fht>;
#    my $last_line = $all_lines[-1];
#    $last_line =~ /^\s*(\d+)/;
#    my $max_burn_in_gen = int($burn_in_fraction * $1); 
#    my @keep_lines = ();
#    my $Tindex_maxLnISW = {};
#    while (@all_lines) {
#       my $line = pop @all_lines; # take the last line
#       my @cols = split(" ", $line);
#       my ($generation, $i_run, $j_T, $k_w, $lnISW) = @cols[0,5,6,7,8];
#       # skip data for generations which are already stored:
#       if ( $generation > $max_burn_in_gen ) {
#          if (exists $Tindex_maxLnISW->{$j_T}) {
#             $Tindex_maxLnISW->{$j_T} = max($Tindex_maxLnISW->{$j_T}, $lnISW);
#          } else {
#             $Tindex_maxLnISW->{$j_T} = $lnISW;
#          }
#          unshift @keep_lines, "$line";
#       } else {
#          last;
#       }
#    }

#    for my $line (@keep_lines) {
#       my @cols = split(" ", $line);
#       my ($generation, $i_run, $j_T, $k_w, $lnISW) = @cols[0,5,6,7,8];
#       #  print "lnISW before offset: ", $lnISW, "  offset: ", $Tindex_maxLnISW->{$j_T}, "\n";
#       $lnISW -= $Tindex_maxLnISW->{$j_T};
#       my $ISW = 1000*exp($lnISW);
#       #  print "lnISW ISW:  $lnISW  $ISW\n";
#       my $sample_set_id = $i_run . "_" . ($j_T % $n_temperatures) . "_" . ($k_w % $n_temperatures); # walkers separate
#       #    print "ssID: $n_temperatures  $sample_set_id.\n";
#       $newick = $cols[9];
#       $newick =~ s/;\s*$//;     # remove final semi-colon.
#       $newick =~ s/^\s+//;      # remove init whitespace.
#       $newick =~ s/:[0-9]+[.][0-9]+(e[-+][0-9]{2,3})?(,|[)])/$2/g; # remove branch lengths
#       $newick =~ s/^\s+//;
#       $newick =~ s/;\s*$//;
#       my $mask = ( 1 << $n_taxa ) - 1; # n_taxa 1's
#       my $split_count = {}; # keys are split bit patterns, values: counts
#       my ( $minlabel, $ordered_newick, $tree_labels, $split_bp ) =
#         Mrbayes::order_newick( $newick, $split_count, $mask, 0 );
#       #    print "$newick   $ordered_newick \n" if($newick ne $ordered_newick);
#       #      $orderednewick_count{$ordered_newick}++;
#       #      $newick_count{$newick}++;
#       $newick = $ordered_newick;
#       if ( $newick =~ s/^\s*[(]1,/'(1,(/ ) { # go from form (1,(...),(...)) with trifurcation, to (1,((...),(...))) with only bifurcations
#          $newick .= ")'";
#       } else {
#          die "newick: $newick doesn't start with '(1,' as expected.\n";
#       }
#       if ( !exists $newick_number_map{$newick} ) {
#          $topology_number++;    # counts distinct topologies

#          $newick_number_map{$newick}          = $topology_number; # 1,2,3,...
#          $number_newick_map{$topology_number} = $newick;
#       }
#       push @generations, $generation;
#       my $ntopo_w = $newick_number_map{$newick} . ":$ISW"; # this is 
#       $Topology_chain_data
#         ->store_data_point( $sample_set_id, $generation, $ntopo_w); # $newick_number_map{$newick} );
#       $Topology_count->{$newick}++;

#       my $the_split = join( ";", keys %$split_count ) . ":$ISW";
#       $Splits_chain_data->store_data_point($sample_set_id, $generation, $the_split );
#    }                            # loop over chunk lines.

#    # $self->{newick_number_map}     = \%newick_number_map;
#    # $self->{number_newick_map}     = \%number_newick_map;
#    # $self->{n_distinct_Topologies} = $topology_number;
#    return ($Topology_chain_data, $Splits_chain_data, \%number_newick_map);
# }

# sub tp_analyze {                # 
#    my $Topology_chain_data = shift;
#    my $Splits_chain_data = shift;
#    my $number_newick = shift;
#    my $file_basename = shift;
#    my $n_runs = shift;
#    my $n_temps = shift;
#    my $n_temps_out = shift;
#    my ($v, $d, $f ) = File::Spec->splitpath($file_basename);
#    my $histogram_filename = $f . '.' . 'Topology_histograms';
#    print "opening file $histogram_filename for writing.\n";
#    open my $fhhist, ">", "$histogram_filename" or die "Couldn't open $histogram_filename for writing.\n";
#    #  print $Topology_chain_data->{histograms}->histogram_string('by_bin_weight', 1.0);
#    print $fhhist $Topology_chain_data->{histograms}->histogram_string('by_bin_weight', 1.0);
#    my @topo_L1_distances = $Topology_chain_data->{histograms}->avg_L1_distance();
#    printf $fhhist ( "# L1:   %6.3f %6.3f %6.3f %6.3f %6.3f \n# Linf: %6.3f %6.3f %6.3f %6.3f %6.3f \n\n", @topo_L1_distances[ 0 .. 4 ], @topo_L1_distances[6..10] );

#    print "Before Ln_distances_x.\n";
#    for my $Tindex (0..$n_temps_out-1) {
#       my $compare_string = '';
#       my @grp_strs = ();
#       for my $Rindex (0..$n_runs-1) {
#          my $rt_bit = $Rindex . '_' . $Tindex . '_';
#          my @rtws = ();
#          for my $Windex (0..$n_temps-1) {
#             my $rtw_bit = $rt_bit . $Windex;
#             push @rtws, $rtw_bit;
#          }
#          push @grp_strs, join(',', @rtws);
#       }
#       $compare_string = join(';', @grp_strs);
#       #print "compare string: $compare_string\n";
#       my ($d1_avg, $d2_avg, $dinf_avg) = $Topology_chain_data->{histograms}->Ln_distances_x($compare_string, 2);
#       print "$Tindex   $d1_avg  $d2_avg  $dinf_avg\n";
#    }



#    # my $number_newick = $self->{number_newick_map};
#    for my $topo_number (1..10000) {
#       if (exists $number_newick->{$topo_number}) {
#          my $newick = $number_newick->{$topo_number};
#          printf $fhhist ("# %4i  %s \n", $topo_number, $newick);
#       } else {
#          last;
#       }
#    }
#    close $fhhist;

#    $histogram_filename = $f . "." . "Splits_histograms";
#    open $fhhist, ">", "$histogram_filename";
#    #  print $fhhist "# After $ngen generations. \n",
#    print $fhhist  $Splits_chain_data->{histograms}->histogram_string('by_bin_weight', 1.0);
#    my @splits_L1_distances = $Splits_chain_data->{histograms}->avg_L1_distance();

#    printf $fhhist ("# L1:   min: %5.3f; q1: %5.3f; median: %5.3f; q3: %6.3f; max %5.3f \n",  @splits_L1_distances[ 0 .. 4 ]);
#    printf $fhhist ("# Linf: min: %5.3f; q1: %5.3f; median: %5.3f; q3: %5.3f; max: %5.3f \n\n", @splits_L1_distances[ 6 .. 10 ] );
#    close $fhhist;

# }

# $self->retrieve_param_samples_from_tp( $tp_filename, $prev_chunk_ngen ); #read in from * .run?.t file

# my %new_param_names = ('Alpha' => 1, 'Pinvar' => 1, 'lnPPr' => 1, 'TreeLength' => 1); #, 'Topology' => 1); 
# my @prams = sort keys %new_param_names;
# for my $the_param (@prams) {
#    #   next if(exists $new_param_names{$the_param});
#    #   print "THE PARAM: $the_param.\n";
#    my $chain_data_obj = $self->{chain_data}->{$the_param};
#    my $histogram_filename =
#      $self->{file_basename} . "." . $chain_data_obj->{parameter_name} . "_histograms";
#    open my $fhhist, ">", "$histogram_filename";
#    print $fhhist $chain_data_obj->{histograms}->histogram_string('by_bin_number');
#    my ($min_L1, $q1_L1, $median_L1, $q3_L1, $max_L1, $max_avg_intercluster_L1 ) = 
#      $chain_data_obj->{histograms}->avg_L1_distance();
#    printf $fhhist ( "# L1 dists: min: %5.3f; q1: %5.3f; median: %5.3f; q3: %5.3f; max: %5.3f; max avg intercluster: %5.3f\n", 
#                     $min_L1, $q1_L1, $median_L1, $q3_L1, $max_L1, $max_avg_intercluster_L1);
#    print $fhhist "# Max Kolmogorov-Smirnov D for parameter ",
#      $chain_data_obj->{parameter_name},
#        " is: ", $chain_data_obj->{histograms}->binned_max_ksd(), "\n\n";
#    close $fhhist;
# }

#    $self->retrieve_topology_samples_from_tp( $tp_filename, $prev_chunk_ngen ); #read in from * .run?.t file
#    my $histogram_filename = $self->{file_basename} . "." . "Topology_histograms";
#    open my $fhhist, ">", "$histogram_filename";
#    print $fhhist "# After $ngen generations. \n",
#      $self->{Topo_chain_data}->{histograms}->histogram_string('by_bin_weight');
#    my @topo_L1_distances = $self->{Topo_chain_data}->{histograms}->avg_L1_distance();
#    printf $fhhist ( "# L1:   %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f \n# Linf: %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f  \n\n", @topo_L1_distances[ 0 .. 11 ] );
#    my $number_newick = $self->{number_newick_map};
#    for my $topo_number (1..10000) {
#       if (exists $number_newick->{$topo_number}) {
#          my $newick = $number_newick->{$topo_number};
#          printf $fhhist ("# %4i  %s \n", $topo_number, $newick);
#       } else {
#          last;
#       }
#    }
#    close $fhhist;

#    $histogram_filename = $self->{file_basename} . "." . "Splits_histograms";
#    open $fhhist, ">", "$histogram_filename";
#    print $fhhist "# After $ngen generations. \n",
#      $self->{Splits_chain_data}->{histograms}->histogram_string('by_bin_weight');
#    my @splits_L1_distances = $self->{Splits_chain_data}->{histograms}->avg_L1_distance();

#    printf $fhhist ("# L1:   min: %5.3f; q1: %5.3f; median: %5.3f; q3: %6.3f; max %5.3f; max avg intercluster: %5.3f \n",  @splits_L1_distances[ 0 .. 5 ]);
#    printf $fhhist ("# Linf: min: %5.3f; q1: %5.3f; median: %5.3f; q3: %5.3f; max: %5.3f; max avg intercluster: %5.3f \n\n", @splits_L1_distances[ 6 .. 11 ] );
#    close $fhhist;
#    return (\@topo_L1_distances, \@splits_L1_distances);
# }

sub analyze_chain_splits{
   my $chain_splits = shift;
   my $n_runs = shift;
   my $n_Ts = shift;
   my $outstring = '';
#print STDERR "Top of analyze_chain_splits. $n_runs  $n_Ts \n";
#print STDERR "rtw:  ", join(" ", keys %$chain_splits), "\n";
   for my $i_run (0.. $n_runs-1) { #
      for my $i_W (0..$n_Ts-1) { #
         my $rtw = $i_run . '_0_' . $i_W; # T index is 0 for now.
#print "RTW: $rtw \n";
         my $gen_splits = $chain_splits->{$rtw}; # this should give a gen:splits hashref
         my  @sgens = sort {$a<=>$b} keys %$gen_splits;
 #        print "# ", join(" ", @sgens[0..10]), "\n";
         my $prev_splits = $gen_splits->{$sgens[0]};
      #   print $prev_splits, "\n"; exit;
         my $prev_gen = 0;
         my $generations = $chain_splits->{generations}->[0];
      #   my @rtw_gens = sort {$a <=> $b} keys
         for my $gen (@sgens) {
        #    print "gen: $gen\n";
            if (exists $gen_splits->{$gen}) {
               my $splits = $gen_splits->{$gen};
               #  if(defined $prev_splits[$i_W]){
               my $TRF = TRF_distance($splits, $prev_splits);
               $outstring .=  "$rtw  " . ($gen - $prev_gen) . "  $TRF\n";
               #   }
               $prev_splits = $splits;
               $prev_gen = $gen;
            }
         }                      # loop over generations
      }                         # loop over walkers
   }                            # loop over runs
   return $outstring;
}

