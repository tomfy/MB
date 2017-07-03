package Histograms;
use strict;
use List::Util qw ( min max sum );
use Statistics::Descriptive;
use lib '/home/tomfy/MB/lib/';
# use Cluster;
# for histogramming Mrbayes chain output
# Histograms can hold several histograms (1 per run)
# and provide some statistics on them, particularly
# how similar or different they are: avg L1 distance,
# Kolmogorov-Smirnov D statistic, etc.
sub new {
   my $class             = shift;
   my $arg               = shift;
   my $default_arguments = {
                            'title'          => 'untitled',
                            'min_gen'        => 0,
                            'n_histograms' => 1,
                            'set_size' => 1
                           };
   my $self = bless $default_arguments, $class;

   foreach my $option ( keys %$arg ) {
      warn "Unknown option: $option in Histograms constructor.\n"
        if ( !exists $self->{$option} );
      if ( defined $arg->{$option} ) { # if arg is undef, leaves default in effect
         $self->{$option} = $arg->{$option};
      }
   }

   my %histograms = ();

   $self->{histograms} = \%histograms; # keys are dataset ids; values are bin:count hashrefs representing histograms
   $self->{sum_histogram} = {};
   $self->{rearr_histograms} =
     {
     }; # same info as histograms, but hashref with keys: bins; values hash refs of setid:counts pairs.
   return $self;
}

sub get_histograms{
   my $self = shift;
   return $self->{histograms}; 
}

sub get_sum_histogram{
   my $self = shift;
   return $self->{sum_histogram};
}

sub get_total_counts {
   my $self = shift;
   return $self->{total_in_histograms};
}

sub get_histogram_counts {
   my $self = shift;
   my $hists = $self->{histograms};
   my $string = '';
   my $total = 0;
   while (my($hid, $cc) = each %$hists) {
      my $hist_sum = sum(values %$cc);
      $total += $hist_sum;
      $string .= "$hid: $hist_sum;  " . join(";", values %$cc) . "\n";
   }
   
   $string .= "xtotal: " . $self->{total_in_histograms} . "  $total";
   return $string;
}

sub set_min_max_bins {
   my $self = shift;
   my ( $min_bin, $max_bin ) = @_;
   $self->{min_bin} = $min_bin;
   $self->{max_bin} = $max_bin;
}
sub adjust_cat_count {
   my $self        = shift;
   my $hist_id = shift;         # 0,1,2,...
   my $category    = shift;
   my $increment   = shift;
   $increment = 1 if ( !defined $increment );
   #   print "in adjust cat count. $hist_id $category $increment. new value:  ";
   $self->{histograms}->{$hist_id} = {} if(!exists $self->{histograms}->{$hist_id});
   $self->{histograms}->{$hist_id}->{$category} += $increment;
   $self->{sum_histogram}->{$category}              += $increment;
   $self->{total_in_histograms}                     += $increment;
   #   print $self->{histograms}->{$hist_id}->{$category}, "\n";
}

sub adjust_val_count {
   my $self = shift;
   $self->adjust_cat_count(@_);
}

sub n_histograms {
   my $self = shift;
   return scalar keys %{ $self->{histograms} };
}

sub title {
   my $self  = shift;
   my $title = shift;
   if ( defined $title ) {
      $self->{title} = $title;
   }
   return $self->{title};
}

sub histogram_string {
   # argument: hashref, keys: category labels; values: array ref holding weights for the histograms
   # (e.g. for different mcmc runs)
   # Output: 1st col: category labels, next n_histograms cols: weights, last col sum of weights over histograms.
   # supply ref to array of categories in desired order, or will sort bins by total weight.
   my $self = shift;
   my $sorting            = shift || 'by_bin_number';
   my $normalization = shift || undef; # undef -> output unnormalized, otherwise normalize s.t. sum = $normalization
   my $max_lines_to_print = shift || 1000; #
   my $histogram_cat_weight = $self->{histograms};
   my $sum_cat_weight       = $self->{sum_histogram};
   my $min_bins_to_show = 50;
   my $min_count_to_show = 0;   #scalar @$histogram_cat_weight;
   my $tail_prob = 0.005; # lump together all least-prob. categories with this total prob.

   my $extra_bit         = '';
   my @sorted_categories = keys %$sum_cat_weight;
   my @sorted_histids = sort keys %{$histogram_cat_weight};
   my %histid_weightsum = map(($_ => 0), @sorted_histids);
   #   print "QWERTYU: ", $self->{title}, "  ", join(";", @sorted_histids), "\n";
   my $total_weight = sum( values %$sum_cat_weight );
   if ( $sorting eq 'by_bin_number' ) {
      @sorted_categories =
        sort { $a <=> $b } @sorted_categories; # keys %$sum_cat_weight;
   } else {                     # sort by avg bin count.
      @sorted_categories =
        sort { 
           if ($sum_cat_weight->{$b} ==  $sum_cat_weight->{$a}) {
              return ($a <=> $b);
           } else {
              return ($sum_cat_weight->{$b} <=> $sum_cat_weight->{$a});
           }
        }
          @sorted_categories;   # keys %$sum_cat_weight;
      while ( scalar @sorted_categories > $min_bins_to_show ) { # remove low-count categories, but leave at least $min_bins_to_show categories
         my $the_cat = $sorted_categories[-1];
         if ( exists $sum_cat_weight->{ $sorted_categories[-1] } ) {
            my $last_bin_count =
              $sum_cat_weight->{ $sorted_categories[-1] };
            if ( defined $last_bin_count ) {
            }

            last
              if ( $last_bin_count >= $min_count_to_show ); # get rid of categories with count of zero.
            pop @sorted_categories;
         } else {
            die "category not found: ", $sorted_categories[-1], "\n";
         }
      }
   }
   for my $histid (@sorted_histids) {
      for my $cat (@sorted_categories) {
         $histid_weightsum{$histid} += $histogram_cat_weight->{$histid}->{$cat} // 0;
      }
   }

   my $total_categories = scalar @sorted_categories;
   my $number_of_bins = ($self->{n_bins})? $self->{n_bins} : $total_categories;
   my $string = "# " . $self->title() . " histogram. ";
   if (defined $self->{n_bins}) {
      $string .= " n_bins: " . $number_of_bins. "  ";
      $string .= "binning_lo_val: " . $self->{binning_lo_val} . "  ";
      $string .= "bin_width: " . $self->{bin_width} . ".\n";
   } else {
      $string .= " categories: $total_categories.\n";
   }
   $string .= (defined $self->{bin_width})? "# bin low edge " : "# bin id   ";
   $string .= join(" ", map(sprintf(" %8s", $_),  @sorted_histids)) . "     sum\n";
   $string .= "#-----------------------------------------------------\n";

   my ( $total_so_far, $shown_total ) = ( 0, 0 );
   my $total_before_this_cat = 0;
   my $count = 0;
   my $sum_total = sum(values %$sum_cat_weight);
   my %hist_tailweight = ();
   my $sum_tailweight = 0;
   for my $cat (@sorted_categories) {
      my $sum_weight = $sum_cat_weight->{$cat};
      $total_before_this_cat = $total_so_far;
      $total_so_far += $sum_weight;

      if ($total_before_this_cat/$sum_total <= (1 - $tail_prob)) {
         my $line_string = '';
         my $cat_string  = $self->category_string($cat);
         $line_string .= $cat_string;

         #   for my $c_w (values %$histogram_cat_weight) {
         for my $histid (@sorted_histids) {
            my $c_w = $histogram_cat_weight->{$histid};
            my $the_weight = (exists $c_w->{$cat})? $c_w->{$cat} : 0;
            #  $histid_weightsum{$histid} += $the_weight;
            if (defined $normalization) {
               $the_weight *= $normalization/$histid_weightsum{$histid};
            }
            $line_string .= sprintf( "%9.4g ", $the_weight );
         }

         $line_string .= sprintf( "%12.7g   %9.4g\n", $sum_weight, $total_so_far/$sum_total );
         $histid_weightsum{sum} += $sum_weight;
         $count++;

         if ( $count > $max_lines_to_print ) {
            $extra_bit = "(Top $max_lines_to_print shown.)";
            last;
         } else {
            $string .= $line_string;
            $shown_total += $sum_weight;
         }
      } else {                  # tail
         for my $histid (@sorted_histids) {
            my $c_w = $histogram_cat_weight->{$histid};
            my $the_weight = (exists $c_w->{$cat})? $c_w->{$cat} : 0;
            if (defined $normalization) {
               $the_weight *= $normalization/$histid_weightsum{$histid};
            }
            $hist_tailweight{$histid} += $the_weight;
           
         }
         #   print STDERR "$cat   $sum_weight $sum_tailweight \n";
         $sum_tailweight += $sum_weight;
         $histid_weightsum{sum} += $sum_weight;
      }
   }                            # loop over categories
   $string .= $self->category_string('other');
   for my $histid (@sorted_histids) {
      $string .=  sprintf( "%9.4g ", $hist_tailweight{$histid});
   }
   $string .= sprintf("%12.7g   %9.4g\n", $sum_tailweight, $total_so_far/$sum_total);
  
   $string .= '# sums:    ';
   for my $histid (@sorted_histids) {
      $string .= sprintf( "%9.4g ", $histid_weightsum{$histid} );
   }
   $string .= sprintf( "%12.7g \n", $histid_weightsum{sum});
   $string .= sprintf("#  total  %10.0g  (%10.0g shown).\n", $total_so_far, $shown_total);
   $string .= "#  in " . $total_categories . " categories. $extra_bit \n";
   return $string;
}

sub category_string {
   my $self = shift;
   my $cat  = shift;
   return sprintf("%6s     ", $cat);
}

sub rearrange_histograms
  { # go from an array of bin:count hashes to hash of bin:(count array).
     my $self = shift;
     my $histograms =
       $self->{histograms};     # hash ref of bin:count hashrefs.
     my $n_histograms     = scalar keys %$histograms;
     my $rearr_histograms = {};

     for my $i_hist ( keys %$histograms ) {
        while (my($bin, $count) = each %{$self->{sum_histogram}}) {
           $rearr_histograms->{$bin}->{$i_hist} = 0;
        }
        my $histogram = $histograms->{$i_hist};
        while ( my ( $bin, $count ) = each %$histogram ) {
           if ( !exists $rearr_histograms->{$bin} ) {
              $rearr_histograms->{$bin} = {};
           }
           $rearr_histograms->{$bin}->{$i_hist} = $histogram->{$bin};
        }
     }
     $self->{rearr_histograms} = $rearr_histograms;
  }

sub _total_variation_distance{
   my $self = shift;
   my $bins = shift;            # array ref to bins in sum_histograms
   my $histogram1 = shift;      # hashref. bin-count pairs
   my $histogram2 = shift;      # hashref. bin-count pairs
   my $sum1 = sum(values %$histogram1);
   return undef if($sum1 <= 0);
   my $sum2 = sum(values %$histogram2);
   return undef if($sum2 <= 0);

   my $sum_abs_diff = 0;
   for my $bin (@$bins) {
      my $count1 = $histogram1->{$bin} || 0;
      my $count2 = $histogram2->{$bin} || 0;
      $sum_abs_diff += abs($count1/$sum1 - $count2/$sum2);
   }
   #  warn "_total_variation_distance. sum1 sum2 should be same in _total_variation_distance: $sum1 $sum2 \n" if($sum1 != $sum2);
   return 0.5*$sum_abs_diff;
}

sub _Ln_distance{
   my $self = shift;
   my $bins = shift;            # array ref to bins in sum_histograms
   my $histogram1 = shift;      # hashref. bin-count pairs
   my $histogram2 = shift;      # hashref. bin-count pairs
   my $sum1 = sum(values %$histogram1);
   return undef if($sum1 <= 0);
   my $sum2 = sum(values %$histogram2);
   return undef if($sum2 <= 0);
   my $n = shift || 2;
   #   my $n_set = $self->{set_size};
   my $n_histograms = scalar keys %{$self->{histograms}};
   #  my $n_sample = $self->get_total_counts() / ($n_set * $n_histograms);
   my $sum_abs_freq_diff_to_n = 0;
   for my $bin (@$bins) {
      my $count1 = $histogram1->{$bin} || 0;
      my $count2 = $histogram2->{$bin} || 0;
      #  $sum1 += $count1; $sum2 += $count2;
      $sum_abs_freq_diff_to_n += (abs($count1/$sum1 - $count2/$sum2))**$n;
   }
   #   warn "_Ln_distance. sum1 sum2 should be same in Ln_distance: $sum1 $sum2 $n_sample \n" if($sum1 != $n_sample*$n_set  or $sum2 != $n_sample*$n_set);

   my $result = (0.5*$sum_abs_freq_diff_to_n)**(1.0/$n);
   return $result;
}

sub _Linfinity_distance{
   my $self = shift;
   my $bins = shift;            # array ref to bins in sum_histograms
   my $histogram1 = shift;      # hashref. bin-count pairs
   my $histogram2 = shift;      # hashref. bin-count pairs
   my $sum1 = sum(values %$histogram1);
   return undef if($sum1 <= 0);
   my $sum2 = sum(values %$histogram2);
   return undef if($sum2 <= 0);
   #  my $n_set = $self->{set_size};
   my $n_histograms = scalar keys %{$self->{histograms}};
   #  my $n_sample = $self->get_total_counts() / ($n_set * $n_histograms);
   my $max_abs_weight_diff = 0;
   for my $bin (@$bins) {
      my $count1 = $histogram1->{$bin} || 0;
      my $count2 = $histogram2->{$bin} || 0;
      #   $sum1 += $count1; $sum2 += $count2;
      $max_abs_weight_diff = max( abs($count1/$sum1 - $count2/$sum2), $max_abs_weight_diff );
   }
   #  warn "_Linfinity_distance. sum1 sum2 should be same in Ln_distance: $sum1 $sum2 $n_sample \n" if($sum1 != $n_sample*$n_set  or $sum2 != $n_sample*$n_set);

   my $result = 0.5*$max_abs_weight_diff;
   return $result;
}

sub tv_distances{
   my $self = shift;
   my $histograms = $self->{histograms};
   my @setids = keys %{$self->{histograms}}; 
   my @labels = sort {$a <=> $b} keys %{$self->{sum_histogram}}; # bin labels
   my @tvds = ();
   my %ij_dist = ();
   my $n_histograms = scalar keys %$histograms;
   die "n_histograms: $n_histograms " . scalar @setids . " not equal.\n" if($n_histograms != scalar @setids);
   for (my $i = 0; $i < $n_histograms; $i++) {
      for (my $j = $i+1; $j < $n_histograms; $j++) {
         my ($setid1, $setid2) = ($setids[$i], $setids[$j]);
         #  print "ZZZ: $i  $j    $setid1  $setid2 \n";
         my $tv_dist = $self->_total_variation_distance(\@labels, $histograms->{$setid1}, $histograms->{$setid2});
         push @tvds, $tv_dist if(defined $tv_dist); # if one of the histograms is empty, tvd is undefined; don't include in @tvds
         $ij_dist{"$i,$j"} = $tv_dist;
      }
   }
   return (\@tvds, \%ij_dist);
}

sub quartiles{   # get the 5 quartiles (i.e. min, q1, median, q3, max)
   my $data = shift;            # ref to array of data values
   my $DStat_obj = Statistics::Descriptive::Full->new();
   $DStat_obj->add_data(@$data);
   my ($min, $q1, $median, $q3, $max) = map( $DStat_obj->quantile($_), (0..4));
   return ($min, $q1, $median, $q3, $max);
}

sub max_avg_intercluster_distance{
   my $self = shift;
   my $n_hists = scalar keys %{$self->{histograms}};
   my @labels = (0..$n_hists-1);
   my $ij_dist = shift;
   #   my $cluster_obj = Cluster->new(\@labels, $ij_dist);
   my ($max_avg_intercluster_dist, $str) = (undef, undef); # $cluster_obj->best_n_partition(); # default: bipartition, max avg inter-cluster dist
   return $max_avg_intercluster_dist;
}

sub dist_stats{
   my $self = shift;
   my $ij_dist = shift;
   my $n_histograms = scalar keys %{$self->{histograms}};
   my @dists = sort {$a <=> $b} values %$ij_dist;
   my @quartiles = quartiles(\@dists);

   my $max_avg_interclstr_dist =  undef; # $self->max_avg_intercluster_distance($ij_dist);

   return (@quartiles, $max_avg_interclstr_dist, $dists[-($n_histograms-1)]);
}

sub tvd_stats{ # statistics summarizing the pairwise total variation distances between chains
   my $self = shift;
   my ($tvds, $ij_dist) = $self->tv_distances();
   # for(keys %$ij_dist){
   #   print "i j tvd: $_  ", $ij_dist->{$_}, "\n";
   # }
   my @tvd_quartiles = quartiles($tvds);
   my $max_avg_interclstr_dist =  $self->max_avg_intercluster_distance($ij_dist);
   return (@tvd_quartiles, $max_avg_interclstr_dist); # min, q1, median, q3, max, intercluster_dist
}

sub maxbindiff_stats{
   my $self = shift;
   my $ij_diff = $self->max_bin_diffs();
   my @mbd_quartiles = quartiles([values %$ij_diff]);
   my $max_avg_interclstr_dist =  $self->max_avg_intercluster_distance($ij_diff);
   return (@mbd_quartiles, $max_avg_interclstr_dist);
}

sub nbadbin_stats{
   my $self = shift;
   my $stringency = shift || 8; # number of bad splits to allow is (N-3)/$stringency
   my $ij_diff = $self->n_bad_bin_diffs();
   my @nbbs = values %$ij_diff;
   my $max_allowed_badbins = $self->{set_size} / $stringency;
   my @nbb_quartiles = quartiles(\@nbbs);
   my $max_avg_interclstr_dist =  $self->max_avg_intercluster_distance($ij_diff);
   my $bad_pairs = 0;
   for (@nbbs) {
      $bad_pairs++ if($_ > $max_allowed_badbins);
   }
   return (@nbb_quartiles, $max_avg_interclstr_dist, $bad_pairs);
}

sub tvd_statsx{
   my $self = shift;
   my $sum_tvd = 0;
   my $sum_tvdsqrd = 0;
   my $max_tvd = -1;
   my $histograms = $self->{histograms};
   my @bins = sort {$a <=> $b} keys %{$self->{sum_histogram}};
   my @tvds = ();
   my @labels = ();
   my %ij_dist = ();
   my $n_histograms = scalar keys %$histograms;
   for (my $i = 0; $i < $n_histograms; $i++) {
      push @labels, $i;
      for (my $j = $i+1; $j < $n_histograms; $j++) {
         my $tv_dist = $self->_total_variation_distance(\@bins, $histograms->{$i}, $histograms->{$j});
         push @tvds, $tv_dist if(defined $tv_dist);
         $sum_tvd += $tv_dist;
         $sum_tvdsqrd += $tv_dist**2;
         $max_tvd = max($max_tvd, $tv_dist);
         $ij_dist{"$i,$j"} = $tv_dist;
      }
   }
   #   my $cluster_obj = Cluster->new(\@labels, \%ij_dist);

   my ($max_avg_intercluster_dist, $str) = (undef, undef); # $cluster_obj->best_n_partition(); # default: bipartition, max avg inter-cluster dist

   my $avg_tvd = $sum_tvd /( $n_histograms*($n_histograms-1)/2 );
   my $rms_tvd = $sum_tvdsqrd /( $n_histograms*($n_histograms-1)/2 )**0.5;
   @tvds = sort {$a <=> $b} @tvds;
   my $x = 3/4 * (scalar @tvds);
   my ($xlo, $xhi) = (int($x), int($x+1));
   my $upper_quartile = $tvds[int($xlo)] * ($xhi - $x) + $tvds[$xhi] * ($x - $xlo);
   my $DStat_obj = Statistics::Descriptive::Full->new() ; # @tvds);
   $DStat_obj->add_data(@tvds);
   return ($avg_tvd, $rms_tvd, $max_tvd, $upper_quartile, $max_avg_intercluster_dist);
}

sub _Ln_distance_old{
   my $self = shift;
   my $bins = shift;
   my $histogram1 = shift;      # hashref. bin-count pairs
   my $histogram2 = shift;      # hashref. bin-count pairs
   my $n = shift || 2;          # L2 distance by default
   my $sum = 0.0;
   #my @sbins = sort {$a <=> $b} @$bins;
   for my $bin (@$bins) {
      my $count1 = $histogram1->{$bin};
      $count1 = 0 if(! defined $count1);
      my $count2 = $histogram2->{$bin};
      $count2 = 0 if(! defined $count2);
      die "Histogram 2 has undefined counts for bin $bin. \n" if( ! defined $count2);
      my $abs_diff = abs($count1 - $count2);
      $sum += $abs_diff**$n;
   }
   return $sum**(1/$n);
}

sub Ln_distances{
   my $self = shift;
   my $n = shift || 2;
   my $avg = 0;
   my $max = -1;
   my $histograms = $self->{histograms};
   my @bins = sort {$a <=> $b} keys %{$self->{sum_histogram}};
   my @distances = ();
   my %ij_dist = ();
   my $n_histograms = scalar keys %$histograms;
   my @histogram_array = keys %$histograms;
   for (my $i = 0; $i < $n_histograms; $i++) {
      for (my $j = $i+1; $j < $n_histograms; $j++) {
         my $Ln_dist = ($n eq 'infinity')?
           $self->_Linfinity_distance(\@bins, $histograms->{$histogram_array[$i]}, $histograms->{$histogram_array[$j]}) :
             $self->_Ln_distance(\@bins, $histograms->{$histogram_array[$i]}, $histograms->{$histogram_array[$j]}, $n);
         push @distances, $Ln_dist if(defined $Ln_dist);
         $ij_dist{"$i,$j"} = $Ln_dist;
         $avg += $Ln_dist;
         $max = max($max, $Ln_dist);
      }
   }
   return \%ij_dist;            # @distances;
}

sub Ln_distances_x{
   my $self = shift;
   my $select_string = shift; # string like: '0_0_1:0.1,0_0_2:0.2,0_0_2:0.13; 0_1_1:0.1,0_1_2:0.2,0_1_2:0.13; ...'
   # sum the comma separated hists (number following : is a weight, which is 1 if absent); semicolon separates groups of histograms, the 
   # distances between the groups (i.e. between the weighted sum of the histograms in the groups) are to be calculated.
   my @categories = keys %{$self->get_sum_histogram()}; # categories (bins) of sum histogram.

   die if (! $select_string);
   my @groups = split(';', $select_string); # each group is a set of one or more histograms to be summed
   #  print "groups:  ", join("   ", @groups), "\n";
   my %group_sumhist = ();
   for my $agroup (@groups) {   # loop over groups to be compared
      my $ghist = {};
      my @subhists = split(',', $agroup); # histograms in the group.
      for my $asubhist (@subhists) { # loop over histograms in group to be added
         my ($histid, $weight) = ($asubhist =~ /(.*)[:](.*)/)? ($1, $2) : ($asubhist, 1.0); # default weight is 1
         for my $cat (@categories) { # loop over all categories
            $ghist->{$cat} += $weight * ($self->get_histograms()->{$histid}->{$cat} // 0);
         }
      }
      $group_sumhist{$agroup} = $ghist;
      # print "Group: $agroup.\n";
      # while(my ($c, $v) = each %$ghist){
      #    print "  $c  $v\n";
      # }
   }  # done getting sum histogram for each group.

   #  print "Categories:  ", join(",", @categories), "\n";
   my ($d1_avg, $d2_avg, $dinf_avg, $count) = (0, 0, 0, 0);
   for (my $i = 0; $i < (scalar @groups - 1); $i++) {
      my $g1 = $groups[$i];
      for (my $j = $i+1; $j < scalar @groups; $j++) {
         my $g2 = $groups[$j];
         my ($d1, $d2, $dinf) = Histograms::_l12inf_distance( \@categories, $group_sumhist{$g1}, $group_sumhist{$g2});
         $d1_avg += $d1; $d2_avg += $d2, $dinf_avg += $dinf; $count++;
         #   print "g1:  $g1;  g2: $g2; d1,2,inf: $d1 $d2 $dinf\n";
      }
   }
   $d1_avg /= $count;
   $d2_avg /= $count;
   $dinf_avg /= $count;
   return ($d1_avg, $d2_avg, $dinf_avg);
}


sub Ln_distances_y{          # compare each group to average of groups
   my $self = shift;
   my $select_string = shift; # string like: '0_0_1:0.1,0_0_2:0.2,0_0_2:0.13; 0_1_1:0.1,0_1_2:0.2,0_1_2:0.13; ...'
   # sum the comma separated hists (number following : is a weight); semicolon separates groups of histograms, the 
   # distances between the groups (i.e. between the weighted sum of the histograms in the groups) are to be calculated.
   my @categories = keys %{$self->get_sum_histogram()}; # categories (bins) of sum histogram.

   die if (! $select_string);
   my @groups = split(';', $select_string); # each group is a set of one or more histograms to be summed
   #  print "groups:  ", join("   ", @groups), "\n";
   my %group_sumhist = ();
   my %overall_sumhist = ();
   for my $agroup (@groups) {   # loop over groups to be compared
      my $ghist = {};
      my @subhists = split(',', $agroup); # histograms in the group.
      for my $asubhist (@subhists) { # loop over histograms in group to be added
         my ($histid, $weight) = ($asubhist =~ /(.*)[:](.*)/)? ($1, $2) : ($asubhist, 1.0); # default weight is 1
         for my $cat (@categories) { # loop over all categories
            my $bin_increment = $weight * ($self->get_histograms()->{$histid}->{$cat} // 0);
            $ghist->{$cat} += $bin_increment;
            $overall_sumhist{$cat} += $bin_increment;
         }
      }
      $group_sumhist{$agroup} = $ghist;
      # print "Group: $agroup.\n";
      # while(my ($c, $v) = each %$ghist){
      #    print "  $c  $v\n";
      # }
   }  # done getting sum histogram for each group.

   #  print "Categories:  ", join(",", @categories), "\n";
   my ($d1_avg, $d2_avg, $dinf_avg, $count) = (0, 0, 0, 0);
   for (my $i = 0; $i < (scalar @groups - 1); $i++) {
      my $g1 = $groups[$i];
      for (my $j = $i+1; $j < scalar @groups; $j++) {
         my $g2 = $groups[$j];
         my ($d1, $d2, $dinf) = Histograms::_l12inf_distance( \@categories, $group_sumhist{$g1}, $group_sumhist{$g2});
         $d1_avg += $d1; $d2_avg += $d2, $dinf_avg += $dinf; $count++;
         #   print "g1:  $g1;  g2: $g2; d1,2,inf: $d1 $d2 $dinf\n";
      }
   }
   $d1_avg /= $count;
   $d2_avg /= $count;
   $dinf_avg /= $count;
   return ($d1_avg, $d2_avg, $dinf_avg);
}

sub max_bin_diffs{ # for each pair of chains, the max over bins of the abs diff in split frequencies
   my $self = shift;
   my $label_weightslist = shift;
   my $set_size = shift || $self->{set_size}; # e.g. if histogramming splits, set_size = N-3,
   if ( !defined $label_weightslist ) {
      $self->rearrange_histograms();
      $label_weightslist = $self
        ->{rearr_histograms}; # hashref. keys are category labels; values are refs to hashes of setid:weight pairs
   }

   # because each topology has N-3 non-terminal splits, all distinct.
   # in the case of histogramming splits, each topology gives N-3 distinct (non-terminal) splits 
   # so e.g. if there are n_chain mcmc chains and you get n_topo topologies from each chain
   # then total_counts = n_chain * n_topo * set_size, and the max possible counts in each bin
   # for the histogram of one chain is n_topo, or total_counts/( n_chain *set_size)

   my $n_histograms           = scalar keys %{$self->{histograms}};
   my $total_counts      = $self->get_total_counts();
   my $max_in_category = $total_counts / ( $n_histograms * $set_size );

   my %ij_maxbindiff = ();
   for (my $i=0; $i<$n_histograms; $i++) {
      for (my $j=$i+1; $j<$n_histograms; $j++) {
         my $ij = "$i,$j";
         $ij_maxbindiff{$ij} = 0;
      }
   }
   my @labels = keys %$label_weightslist; #
   foreach my $label (@labels) { # loop over categories (bins)
      my $setid_weight = $label_weightslist->{$label};
      next if(sum(values %$setid_weight) <= 0);
      my @setids = sort keys %{$setid_weight};
      die "n_histograms: $n_histograms and number of setids " . scalar @setids, " not equal.\n" if($n_histograms != scalar @setids);
      for (my $i=0; $i<$n_histograms; $i++) {
         for (my $j=$i+1; $j<$n_histograms; $j++) {
            my ($setid1, $setid2) = ($setids[$i], $setids[$j]);
            my $ij = "$setid1,$setid2";
            my $ijdiff = abs($setid_weight->{$setid1} - $setid_weight->{$setid2}) / $max_in_category;
            if ($ijdiff > $ij_maxbindiff{$ij}) {
               $ij_maxbindiff{$ij} = $ijdiff;
            }
            ;
         }
      }
   }
   return \%ij_maxbindiff;
}

sub n_bad_bin_diffs{ # for each pair of chains, the number of bins with frequency diff > threshold
   my $self = shift;
   my $threshold = shift || 0.2;
   my $label_weightslist = shift;
   my $set_size = shift || $self->{set_size}; # e.g. if histogramming splits, set_size = N-3,
   if ( !defined $label_weightslist ) {
      $self->rearrange_histograms();
      $label_weightslist = $self
        ->{rearr_histograms}; # hashref. keys are category labels; values are refs to arrays of weights
   }

   # because each topology has N-3 non-terminal splits, all distinct.
   # in the case of histogramming splits, each topology gives N-3 distinct (non-terminal) splits 
   # so e.g. if there are n_chain mcmc chains and you get n_topo topologies from each chain
   # then total_counts = n_chain * n_topo * set_size, and the max possible counts in each bin
   # for the histogram of one chain is n_topo, or total_counts/( n_chain *set_size)

   my $n_histograms           = scalar keys %{$self->{histograms}};
   my $total_counts      = $self->get_total_counts();
   my $max_in_category = $total_counts / ( $n_histograms * $set_size );

   my %ij_nbadbins = ();
   for (my $i=0; $i<$n_histograms; $i++) {
      for (my $j=$i+1; $j<$n_histograms; $j++) {
         my $ij = "$i,$j";
         $ij_nbadbins{$ij} = 0;
      }
   }
   my @labels = keys %$label_weightslist; #
   foreach my $label (@labels) {          # loop over categories
      #  my @weights = @{ $label_weightslist->{$label} };
      my $setid_weight = $label_weightslist->{$label};
      next if(sum(values %$setid_weight) <= 0);
      my @setids = sort keys %{$setid_weight};
      for (my $i=0; $i<$n_histograms; $i++) {
         for (my $j=$i+1; $j<$n_histograms; $j++) {
            my ($setid1, $setid2) = ($setids[$i], $setids[$j]);
            my $ij = "$setid1,$setid2";
            my $ijdiff = abs($setid_weight->{$setid1} - $setid_weight->{$setid2}) / $max_in_category;
            if ($ijdiff > $threshold) {
               $ij_nbadbins{$ij}++;
            }
         }
      }
   }
   return \%ij_nbadbins;
}

sub get_avg_histogram{
   my $self = shift;
   my $avg_histogram = {};
   my $n_histograms = scalar keys %{$self->{histograms}};
   for my $hist (values %{$self->{histograms}}) {
      my $sum = sum(values %$hist);
      while (my ($bin, $v) = each %$hist) {
         $avg_histogram->{$bin} += ($v // 0)/($sum * $n_histograms);
      }
   }
   $self->{avg_histogram} = $avg_histogram;
   return $avg_histogram;
}

sub get_l12s{            # get l1 and l2 distances from mean histogram
   my $self = shift;
   my @l1s = ();
   my @l2s = ();
   # if(!exists $self->{avg_histogram}){
   #    $self->get_avg_histogram();
   # }elsif(sum(values %{$self->{avg_histogram}}) != 1){
   #  print "sum of avg histogram: ", sum(values %{$self->{avg_histogram}}), "\n";
   $self->get_avg_histogram();
   #  }
   #  print "sum of avg histogram: ", sum(values %{$self->{avg_histogram}}), "\n"; #exit;
   my $n_histograms = scalar keys %{$self->{histograms}};
   my @skeys = sort keys  %{$self->{histograms}};
   for my $histid (@skeys) {
      my $hist = $self->{histograms}->{$histid};
      my $sum = sum(values %$hist);
      my ($l1, $l2) = (0, 0);
      while (my ($bin, $v) = each %$hist) {
         my $abs_diff = abs($v/$sum - $self->{avg_histogram}->{$bin});
         $l1 += $abs_diff;
         $l2 += ($abs_diff)**2;
      }
      $l1 *= 0.5;               # so max possible l1 is 1 (i.e. tvd)
      $l2 = sqrt(0.5*$l2);      # so max possible l2 is 1
      push @l1s, $l1;
      push @l2s, $l2;
   }
   return (\@l1s, \@l2s);
}

sub get_l12s_x{          # get l1 and l2 distances from mean histogram
   my $self = shift;
   my $histograms = shift // $self->{histograms};
   my $sum_histogram = shift;
   my @l1s = ();
   my @l2s = ();
   # if(!exists $self->{avg_histogram}){
   #    $self->get_avg_histogram();
   # }elsif(sum(values %{$self->{avg_histogram}}) != 1){
   #  print "sum of avg histogram: ", sum(values %{$self->{avg_histogram}}), "\n";
   if (!defined $sum_histogram) {
      $self->get_avg_histogram();
      $sum_histogram = $self->{avg_histogram};
   }
      #  print "sum of avg histogram: ", sum(values %{$self->{avg_histogram}}), "\n"; #exit;
      #  my $n_histograms = scalar keys %{$histograms};
      my @skeys = sort keys  %{$histograms};
      for my $histid (@skeys) {
         my $hist = $histograms->{$histid};
         my $sum = sum(values %$hist);
         my ($l1, $l2) = (0, 0);
         while (my ($bin, $v) = each %$hist) {
            my $abs_diff = abs($v/$sum - $sum_histogram->{$bin});
            $l1 += $abs_diff;
            $l2 += ($abs_diff)**2;
         }
         $l1 *= 0.5;            # so max possible l1 is 1 (i.e. tvd)
         $l2 = sqrt(0.5*$l2);   # so max possible l2 is 1
         push @l1s, $l1;
         push @l2s, $l2;
      }
      return (\@l1s, \@l2s);
   }

   sub avg_L1_distance { # just find for histograms as given, no rebinning.
      my $self              = shift;
      my $label_weightslist = shift;
      if ( !defined $label_weightslist ) {
         $self->rearrange_histograms();
         $label_weightslist = $self
           ->{rearr_histograms}; # hashref. keys are category labels; values are refs to arrays of weights
      }
      my $set_size = shift || $self->{set_size}; # e.g. if histogramming splits, set_size = N-3,
      # because each topology has N-3 non-terminal splits, all distinct.
      # in the case of histogramming splits, each topology gives N-3 distinct (non-terminal) splits 
      # so e.g. if there are n_chain mcmc chains and you get a sample of size N from each chain
      # then total_counts = n_chain * N * set_size, and the max possible counts in each bin
      # for the histogram of one chain is N, or total_counts/( n_chain *set_size)

      my $n_histograms           = scalar keys %{$self->{histograms}};

      my $total_counts      = $self->get_total_counts();
      my ($above_threshold_string1, $above_threshold_string2, $above_threshold_string3) = ('', '', '');
      my $max_in_category = $total_counts / ( $n_histograms * $set_size );

      # one weight for each histogram being compared (e.g. 1 for each MCMC chain)
      my @thresholds = (0.1, 0.2, 0.4);
      my %threshold_count1 = (); my %threshold_count2 = (); my %threshold_count3 = ();
      for (@thresholds) {
         $threshold_count1{$_} = 0;
         $threshold_count2{$_} = 0;
         $threshold_count3{$_} = 0;
      }

      #counts categories with abs diff > the keys
      my ($avg_L1_distance, $max_range ) = ( 0, 0, -1 );
      my @labels = keys %$label_weightslist; #
      # my %sumw_bsos = ();
      my ($sum_ranges, $sumsq_ranges, $src) = (0, 0, 0);

      if ( $n_histograms > 1 ) { #
         my $sum_absdiffs = 0;
         my $counts_each_run = $total_counts / $n_histograms;
         my %label_range   = ();
         foreach my $label (@labels) { # loop over categories
            my @weights = sort { $b <=> $a } values %{ $label_weightslist->{$label} };
            next if(sum(@weights) <= 0);

            my ($mean, $variance) = mean_variance(\@weights);
            my $binomial_bin_variance = max(0.0, $mean*(1 - $mean/$max_in_category)); # variance of binomial distribution with mean $mean
         
            my $this_label_range = (max(@weights) - min(@weights)) / $max_in_category;
            $sum_ranges += $this_label_range;
            $sumsq_ranges += $this_label_range**2;
            $src++;
            $label_range{$label} = $this_label_range;

            #   print "variance: $variance  max in category: $max_in_category  mean: $mean . $binomial_bin_variance\n";
            my ($obs_bin_stddev, $binomial_bin_stddev) = (sqrt($variance)/$max_in_category, sqrt($binomial_bin_variance/$max_in_category));
            # $sumw_bsos{sum(@weights)} = [$binomial_bin_stddev, $obs_bin_stddev];
            for ( @thresholds ) {
               $threshold_count1{$_}++  if ($this_label_range > $_ );
            }
            for ( @thresholds ) {
               $threshold_count2{$_}++ if ($obs_bin_stddev > ($_ * $binomial_bin_stddev));  
            }
            for (my $i = 0; $i < $n_histograms; $i++) {
               for (my $j = 0; $j < $n_histograms; $j++) {
                  my $absdiff = abs($weights[$i] - $weights[$j]) / $max_in_category;
                  for (@thresholds) {
                     $threshold_count3{$_}++ if( $absdiff > $_);
                  }
               }
            }
            $max_range = max( $max_range, $this_label_range );
            my $coeff = $n_histograms - 1;
            for my $histogram_weight (@weights) { # loop over runs
               $sum_absdiffs += $coeff * $histogram_weight;
               #	$histogram_weight;	# accumulating counts in all runs
               $coeff -= 2;
            }
         }
         $above_threshold_string1 =
           join( " ", map( $threshold_count1{$_}, @thresholds));
         $above_threshold_string2 =
           join( " ", map( $threshold_count2{$_}, @thresholds ) );
         $above_threshold_string3 =
           join( " ", map( $threshold_count3{$_}, @thresholds ) );
         $avg_L1_distance = $sum_absdiffs / ( $total_counts * ( $n_histograms - 1 ) );
      }                                        # loop over histograms
      my @tvd_sstats = $self->tvd_stats();     # 6 numbers
      my @mbd_sstats = $self->maxbindiff_stats(); # 6 numbers
      my @nbb_sstats = $self->nbadbin_stats();    # 7 numbers
      my ($xtvds, $xij_dist) = $self->tv_distances();
      my $ij_L1d = $self->Ln_distances(1);
      my $ij_L2d = $self->Ln_distances(2);
      my $ij_Linfd = $self->Ln_distances('infinity');
      my $ij_nbbd = $self->n_bad_bin_diffs();
      # my @L1_stats = $self->dist_stats($ij_L1d);
      # my @L2_stats = $self->dist_stats($ij_L2d);
      # my @Linf_stats = $self->dist_stats($ij_Linfd);
      # my @nbbd_stats = $self->dist_stats($ij_nbbd);
      # print "tvdstats:  ", join("; ", @tvd_sstats), "\n";
      # print "L1_stats:  ", join(": ", @L1_stats), "\n";
      # print "L2_stats:  ", join(": ", @L2_stats), "\n";
      # print "mbdstats:  ", join(": ", @mbd_sstats), "\n";
      # print "Linf_stats:  ", join(": ", @Linf_stats), "\n";
      # print "nbb_sstats :  ", join("; ", @nbb_sstats), "\n";
      # print "nbbd_stats:  ", join(": ", @nbbd_stats), "\n";
      #print join("; ", @$xtvds), "\n";
      #print join(": ", value), "\n";
      #print "nbbstats: ", join("; ", @nbb_sstats), "\n";
      my $n_bad_pairs = pop @nbb_sstats; 

      #print"nbbstats: ", join("; ", @nbb_sstats), "  xxx   ", $n_bad_pairs, "\n";
      my $above_threshold_string4 = join(" ", @nbb_sstats);
      #return (\@L1_stats, \@L2_stats, \@Linf_stats, \@nbbd_stats);
      return ( @tvd_sstats, @mbd_sstats, $above_threshold_string1, $above_threshold_string2, $above_threshold_string3, $above_threshold_string4, $n_bad_pairs );
   }

   sub four_distance_stats{
      my $self = shift;
      my $ij_L1d = $self->Ln_distances(1);
      my $ij_L2d = $self->Ln_distances(2);
      my $ij_Linfd = $self->Ln_distances('infinity');
      my $ij_nbbd = $self->n_bad_bin_diffs();
      my @L1_stats = $self->dist_stats($ij_L1d);
      my @L2_stats = $self->dist_stats($ij_L2d);
      my @Linf_stats = $self->dist_stats($ij_Linfd);
      my @nbbd_stats = $self->dist_stats($ij_nbbd);
      return (\@L1_stats, \@L2_stats, \@Linf_stats, \@nbbd_stats);
   }

   sub binned_max_ksd {
      my $self         = shift;
      my $n_histograms = scalar keys %{$self->{histograms}};
      my $max_ksd      = 0;
      my %cume_probs  = ();     #  = ( (0) x $n_histograms );

      my @setids = keys %{$self->{histograms}};
      for (@setids) {
         $cume_probs{$_} = 0;
      }
      my $total_counts = sum( values %{ $self->{histograms}->{$setids[0]}} ); # counts in one representative histogram

      for my $bin ( $self->{min_bin} .. $self->{max_bin} ) {
         my $histograms = $self->{histograms};

         for my $a_setid (keys %$histograms ) {
            my $b_c = $histograms->{$a_setid};
            $cume_probs{$a_setid} += $b_c->{$bin};
         }
         my $cdf_range = max(values %cume_probs) - min(values %cume_probs);
         if ( $cdf_range > $max_ksd ) {
            $max_ksd = $cdf_range;
         }
      }
      return ($total_counts > 0)? $max_ksd/$total_counts : '---';
   }

   sub mean_variance{
      my $x = shift;            # array ref of numbers.
      my ($count, $sum_x, $sum_xsqr) = (scalar @$x, 0, 0);
      return (undef, undef) if($count == 0);
      for (@$x) {
         $sum_x += $_;
         $sum_xsqr += $_*$_;
      }
      my $mean = $sum_x/$count;
      my $variance = ($count > 1)? ($count/($count-1)) * $sum_xsqr/$count - $mean**2 : undef; # 'bessel correction n/(n-1) applied -> unbiased estimator of population variance.
      #  my $stderr = ($count > 1)? sqrt($variance/$count) : undef;
      return ($mean, $variance); #, $stderr);
   }

   # sub minweight_L1 {
   #     my $self = shift;
   #     my $minweight = shift || 0.02;
   #     return $self->avg_L1_distance( $self->minweight_rebin($minweight) )
   #       ;    # $self->minweight_rebin($minweight));
   # }

   # sub minweight_rebin {

   # # input here is already binned
   # # the idea here is to make each bin have at least some fraction of total weight (1% is default)
   # # by possibly lumping together some bins.
   #     my $self = shift;
   #     my $target_bin_weight = shift || 0.01;
   #     $self->rearrange_histograms();
   #     my $label_weightslist = $self->{rearr_histograms};

   # #print "in minweight_rebin. labels: ", join(", ", keys %$label_weightslist), "\n";

   #     my %label_sumweights = ();
   #     while ( my ( $l, $ws ) = each %$label_weightslist ) {
   #         $label_sumweights{$l} = sum(@$ws);
   #     }
   #     my @sorted_labels =
   #       sort {
   #         $label_sumweights{$a} <=> $label_sumweights{$b}
   #       }    # sort by weight; small to large
   #       keys %$label_weightslist;

   #  #print "in minweight_rebin. sorted labels: ", join(", ", @sorted_labels), "\n";

   #     my $total_hits =
   #       sum( map( @{ $label_weightslist->{$_} }, @sorted_labels ) );
   #     my $run0_hits = sum( map( $label_weightslist->{$_}->[0], @sorted_labels ) );
   #     my $n_histograms = scalar @{ $label_weightslist->{ $sorted_labels[0] } };

   #     my $result       = {};
   #     my $cume_weight  = 0;
   #     my @cume_weights = ( (0) x $n_histograms );
   #     my $cume_label   = '';
   #     foreach my $label (@sorted_labels) {    # loop over categories
   #         my @weights = @{ $label_weightslist->{$label} };
   #         @cume_weights = map { $cume_weights[$_] + $weights[$_] } 0 .. $#weights;
   #         my $weight =
   #           sum(@weights); # number of hits for this categories, summed over runs.
   #         $cume_weight += $weight;
   #         $cume_label .= $label . '_';
   #         if ( $cume_weight >= $target_bin_weight * $total_hits ) {
   #             my @copy = @cume_weights;
   #             $cume_label =~ s/_$//;
   #             $result->{$cume_label} = \@copy;
   #             $cume_weight = 0;
   #             @cume_weights = ( (0) x $n_histograms );
   #             $cume_label = '';
   #         }
   #     }
   #     return $result;
   # }

   # sub populate {
   #     my $self                = shift;
   #     my $histogram_gen_value = shift
   #       ; # ref to array of gen/value hashrefs. $rgv->[0]->{120} is value for run 0, gen 120
   #     my $min_gen = shift || 0;

   # #  while (my ($run, $g_v) = each @$histogram_gen_value) { # can  use this if perl 5.12 or later.
   #     for my $run ( 0 .. $#$histogram_gen_value ) {
   #         my $g_v = $histogram_gen_value->{$run};
   #         for my $g ( sort { $a <=> $b } keys %$g_v ) {
   #             my $v = $g_v->{$g};
   #             if ( $g >= $min_gen ) {
   #                 $self->{histograms}->[$run]->{$v}++;
   #                 $self->{sum_histogram}->{$v}++;
   #             }
   #         }
   #     }
   # }

   # sub get_count{
   #   my $self = shift;
   #   my $hist_number = shift;
   #   my $category = shift;
   #   return $self->{histograms}->[$hist_number]->{$category};
   # }


   sub _l12inf_distance{
      my $categories = shift;
      my $hist1 = shift;
      my $hist2 = shift;
      my ($l1dist, $l2dist, $linfdist) = (0, 0, 0);
      my ($sum1, $sum2) = (sum(values %$hist1), sum(values %$hist2));
      return undef if($sum1 == 0  or  $sum2 == 0);
      for my $c (@$categories) {
         my $abs_delta =  abs( ($hist1->{$c} // 0)/$sum1 -  ($hist2->{$c} // 0)/$sum2 );
         $linfdist = max($linfdist, $abs_delta);
         $l1dist += $abs_delta;
         $l2dist += $abs_delta**2;
      }
      # normalized s.t. with hist1 hist2 normalized to sum to 1, these distances will be in interval [0,1], so l1dist is total variation dist.
      $l1dist *= 0.5;
      $l2dist = (0.5*$l2dist)**0.5;
      return ($l1dist, $l2dist, $linfdist);
   }


   package BinnedHistograms;
   use strict;
   use List::Util qw ( min max sum );
   use base qw/ Histograms /;

   sub new {
      my $class = shift;
      my $self  = $class->SUPER::new(@_);
      $self->{n_bins} //= 24;
      $self->{binning_lo_val} //= 0.0;
      $self->{bin_width} //= 0.5;
      # warn "BinnedHistogram constructed with undefined binning parameters.\n"
      #  if(!defined $self->{binning_lo_val} or !defined $self->{bin_width});
      print "in BinnedHistograms constructor: ", $self->{n_bins}, "  ", $self->{binning_lo_val}, "  ", $self->{bin_width}, "\n";
      return $self;
   }

   sub binning_info_string {
      my $self   = shift;
      my $string = '';
      $string .=
        "N bins: " . $self->{n_bins} // 'undef'
          . "Binning low value: " 
            . $self->{binning_lo_val} // 'undef'
              . ". Bin width: "
                . $self->{bin_width} // 'undef' . " \n";
      return $string;
   }

   sub adjust_val_count {
      my $self        = shift;
      my $set_id = shift;       # 0,1,2,...
      my $value       = shift;
      my $increment   = shift;
      $increment = 1 if ( !defined $increment );

      my $category = $self->bin_the_point($value);

      $self->{histograms}->{$set_id}->{$category} += $increment;
      $self->{sum_histogram}->{$category}              += $increment;
      $self->{total_in_histograms}                     += $increment;
   }

   sub bin_the_point {
      my $self  = shift;
      my $value = shift;
      # print "top of bin the point.\n";
      # print "n bins: " , $self->{n_bins} // 'undef', "  ";
      # print "low val: ", $self->{binning_lo_val} // 'undef', "  ";
      # print "bin width: ", $self->{bin_width} // 'undef', "  ";
      die "In bin_the_point. n_bins, binning_lo_val or bin_width not defined. \n"
        if ( !defined $self->{n_bins}
             or !defined $self->{binning_lo_val}
             or !defined $self->{bin_width} );
      my $bin = int( ( $value - $self->{binning_lo_val} ) / $self->{bin_width} );
   
      $bin = max( $bin, 0 );                # bin 0 is 'underflow'
      $bin = min( $bin, $self->{n_bins} - 1 ); # bin n_bins-1 is 'overflow'
      #  print " bin: $bin \n";
      $self->{min_bin} = min( $bin, $self->{min_bin} );
      $self->{max_bin} = max( $bin, $self->{max_bin} );
      return $bin;
   }

   sub category_string {
      my $self    = shift;
      my $cat     = shift;
      my $lo_edge = $self->{binning_lo_val} + $cat * $self->{bin_width};
      return sprintf( "%11.5f     ", $lo_edge );
   }

   # sub populate {
   #     my $self                = shift;
   #     my $histogram_gen_value = shift
   #       ; # ref to array of gen/value hashrefs. $rgv->[0]->{120} is value for run 0, gen 120
   #     my $min_gen = shift || 0;

   # #  while (my ($run, $g_v) = each @$histogram_gen_value) { # can  use this if perl 5.12 or later
   #     for my $run ( 0 .. $#$histogram_gen_value ) {
   #         my $g_v = $histogram_gen_value->[$run];
   #         for my $g ( sort { $a <=> $b } keys %$g_v ) {
   #             my $v = $g_v->{$g};
   #             if ( $g >= $min_gen ) {
   #                 $v = $self->bin_the_point($v);
   #                 $self->{histograms}->[$run]->{$v}++;
   #                 $self->{sum_histogram}->{$v}++;
   #             }
   #         }
   #     }
   # }

   # non-method subroutines





   1;
