package ChainData;
use strict;
use List::Util qw ( min max sum );
use lib '/home/tomfy/MB/lib';
use Histograms;
# for handling Mrbayes parameter output

# param name (topology, logL, treelength, etc.
# generations list, and pre-burn-in gen range and post-burn-in range
# array of gen/value hashrefs, one for each run

#
# max_gen     largest generation number so far
# generations  not used ???
# run_gen_value   $self->{run_gen_value}->[$run]->{$gen} is the value (numerical param, or topology, ...)
# at generation $gen of run $run
# n_runs
# binnable
# histograms
#
#

sub  new {
   my $class = shift;
   my $arg = shift;
   my $default_arguments = {
                            'parameter_name' => 'unnamed parameter',
                            'n_runs' => undef,
                            'n_temperatures' => undef,
                            'n_temperatures_out' => undef,
                            'burn_in_fraction' => 0.1,
                            'gen_spacing' => undef,
                            'binnable' => 1, # by default, assume this is continuous data which should be binned.
                            'next_bin_gen' => 0,
                            'set_size' => 1
                           };
   my $self = bless $default_arguments, $class;

   foreach my $option (keys %$arg) {
      warn "Unknown option: $option in ParamData constructor.\n" if(!exists $self->{$option});
      if (defined $arg->{$option}) { # if arg is undef, leaves default in effect
         $self->{$option} = $arg->{$option};
         #      print STDERR "In ChainData constructor. Setting option $option to ", $arg->{$option}, "  ", $self->{$option},"\n";
      }
   }
   my $n_rtw = $self->{n_runs} * $self->{n_temperatures_out} * $self->{n_temperatures}; # keep track of each run/temperature/walker combination separately
   $self->{max_gen} = 0;
   $self->{generations} = [ ];	# ref to array of n_runs arrayrefs
   $self->{setid_gen_value} = {}; #set id would be something like 0_1_3 for run 0, temperatur 1, walker 3 
   for (0..$self->{n_runs}-1) {
      push @{$self->{generations}}, []; # so $self->{generations} is ref to an array of array refs.
   }

   if ($self->{binnable}) {
      $self->{histograms} = BinnedHistograms->new({title => $self->{parameter_name}, set_size => $self->{set_size}});
   } else {
      $self->{histograms} = Histograms->new({title => $self->{parameter_name}, set_size => $self->{set_size}});
   }
   return $self;
}

sub size_info{
   my $self = shift;
   my $string = "parameter: " . $self->{parameter_name} . "\n";
   my $total_counts = 0;
   while (my($setid, $gv) = each %{$self->{setid_gen_value}}) {
      my @gens = sort {$a<=>$b} keys %$gv;
      $total_counts += scalar @gens;
      $string .= "setid, n points:  $setid  " . scalar @gens . "  " . join(",", @gens) . "\n";
   }
   $string .= "total histogram counts: " . $self->{histograms}->get_histogram_counts() . "  $total_counts.\n";
   return $string;
}


sub store_data_point{		# store a new data point
   # can be single number, or string of ';' separated numbers: '345;277;921'
# optionally followed by : and a weight. e.g.  234;55;9093:0.4556
# which means each of the 3 numbers 234, 55, and 9093 is to be stored with weight 0.4556
# this is used with splits; for each state there is a weight, and a tree topology or equivalently
# a set of splits, each of which gets the same weight.
   my $self = shift;
   my ($setid, $gen, $value) = @_;
   my ($val, $weight) = ($value =~ /(\S+)[:](\S+)/)? ($1, $2) : ($value, 1);
   # if ($val =~ /(\S+)[:](\S+)/) { # get weight if present (else weight is 1)
   #    $val = $1;
   #    $weight = $2;
   # }
   #    print "$setid $gen $val $weight.\n";  
   $self->{max_gen} = max($gen, $self->{max_gen});
   $self->{setid_gen_value}->{$setid} = {} if(!exists $self->{setid_gen_value}->{$setid});
   $self->{setid_gen_value}->{$setid}->{$gen} = $value; # store a value which includes weight
   for my $v (split(";", $val)) {
   #   my ($v,$weight) = ($vw =~ /(\S+)[:](\S+)/)? ($1, $2) : ($vw, 1);
      if (!($self->{binnable})) {
         $self->{histograms}->adjust_val_count($setid, $v, $weight);
      } elsif (defined $self->{histograms}->{binning_lo_val}) {
         $self->{histograms}->adjust_val_count($setid, $v, $weight);
      } else { # else it's binnable, but binning params not set yet; don't adjust val
      }
   }
#   print "new total histogram counts after storing data point:\n", $self->{parameter_name}, "   setid  $setid   ", $self->{histograms}->get_histogram_counts(), "\n";
#   print $self->size_info();
}

sub delete_data_point{
   # if argument is e.g. '1;2;5', splits on ';' and deletes
   # 1, 2 and 5 from histogram. (useful with splits)
   # as with store_data_point, can have weight
   my $self = shift;
   my ($setid, $gen) = @_;
   my $value = $self->{setid_gen_value}->{$setid}->{$gen};
my ($val, $weight) = ($value =~ /(\S+)[:](\S+)/)? ($1, $2) : ($value, 1);
   # my $weight = 1;
   #    if($value =~ /(\S+)[:](\S+)/){
   #       $value = $1;
   #       $weight = $2;
   #    }
#   print "In delete_data_point. setid, gen to delete: $setid, $gen.\n";
   delete $self->{setid_gen_value}->{$setid}->{$gen};
   for my $v (split(";", $val)) {
   #   my ($v,$weight) = ($vw =~ /(\S+)[:](\S+)/)? ($1, $2) : ($vw, 1);
#      print "before adjust_val_count.\n";
      $self->{histograms}->adjust_val_count($setid, $v, -1 * $weight);
   }
#   print "new total histogram counts after deleteing data point:\n", $self->{parameter_name}, "   setid  $setid   \n", $self->{histograms}->get_histogram_counts(), "\n";
#   print $self->size_info();
}

sub delete_low_gens{
   my $self = shift;
   my $max_gen_to_delete = shift;
   my $delta_gen = $self->{gen_spacing};
   # for my $setid (keys %{$self->{setid_gen_value}}) {
   while (my($setid, $gv) = each %{$self->{setid_gen_value}}) {
   #   print "setid: $setid  maxgentodelete: $max_gen_to_delete.  ", join(" ", sort {$a<=>$b} keys %$gv), "\n";
      for (my $g = $max_gen_to_delete; $g >= 0; $g -= $delta_gen) {
       #  print "XYZ:  $g $max_gen_to_delete; ", ;
         if (exists $self->{setid_gen_value}->{$setid}->{$g}) {
            my $value =  $self->{setid_gen_value}->{$setid}->{$g};
          #  print "before delete_data_point.\n";
            $self->delete_data_point($setid, $g);
         } else {
            next;
         }
      }
   #   print "\n";
   }
}

sub delete_pre_burn_in{
   my $self = shift;
   my $max_pre_burn_in_gen = int($self->{max_gen}*$self->{burn_in_fraction});
 #  print "In delete_pre_burn_in; ", $self->{max_gen}, "  ", $self->{burn_in_fraction}, "  max_pre_burn_in_gen: $max_pre_burn_in_gen \n";
 #  print "QQQQQQ:", $self->size_info(), "\n";
   $max_pre_burn_in_gen = $self->{gen_spacing}*int($max_pre_burn_in_gen/$self->{gen_spacing});
   $self->delete_low_gens($max_pre_burn_in_gen);
  # print "RRRRRR: ", $self->size_info(), "\n";
}

sub get_set_data{ # returns hashref with generation/parameter value pairs
   my $self = shift;
   my $setid = shift;		# 0-based
   return $self->{setid_gen_value}->{$setid};
}

sub get_setid_gen_value{ # returns array ref holding hashrefs with generation/parameter value pairs
   my $self = shift;
   return $self->{setid_gen_value};
}

sub get_param_name{
   my $self = shift;
   return $self->{parameter_name};
}


sub get_binning_parameters{
   my $self = shift;
   my $n_bins = shift || 40;
   my $tail_p = shift || 0.05;
   my $tail_d = shift || 0.25;
   # Find the range of data after excluding the lowest $tail_p and the highest $tail_p
   # add $tail_d of this range on each end, divide this new range into $n_bin bins.
   # e.g. 

   #  my @rgv = @{$self->{run_gen_value}};
   my %sgv = %{$self->{setid_gen_value}};

   my @values = ();
   
   for my $g_v (values %sgv) {
      @values = (@values, values %$g_v);
   }
   @values = sort {$a <=> $b} @values;
 #  print "in get binning params: \n"; # , join(" ", @values), "\n";
   my $size = scalar @values;
   my ($lo_val, $hi_val) = ($values[int($tail_p*$size)], $values[int((1 - $tail_p)*$size)]);
   my $central_range = $hi_val - $lo_val;

   $lo_val -= $tail_d*$central_range;
   $hi_val += $tail_d*$central_range;
   if ($lo_val < $values[0]) {
      $lo_val = 0.5*sum($values[0], $lo_val);
   }
   ;
   if ($hi_val > $values[-1]) {
      $hi_val = 0.5*sum($values[-1], $hi_val);
   }
   ;
   if ($lo_val > 0  and   $lo_val < 0.1*$hi_val) {
      $lo_val = 0;
   }
   my $bin_width = ($hi_val - $lo_val)/$n_bins;
   $bin_width = 1 if($bin_width == 0);
   return ($lo_val, $bin_width);
}

sub bin_the_data{
   my $self = shift;
   my $n_bins = shift || 24;
   my $tail_p = shift || 0.05;
   my $tail_d = shift || 0.3;
   my ($lo_val, $bin_width) = $self->get_binning_parameters($n_bins, $tail_p, $tail_d);
 #  print "In bin_the_data. n bins, lo val, bin width:  $n_bins  $lo_val  $bin_width \n";
   my $the_histograms =  $self->{histograms} = 
     BinnedHistograms->new({title => $self->{parameter_name}, n_bins => $n_bins, 
                            binning_lo_val => $lo_val, bin_width => $bin_width});
   my ($min_bin, $max_bin) = ($n_bins, 0);
   my %sgv = %{$self->{setid_gen_value}};
   for my $setid (keys %sgv) {
      my $g_v = $sgv{$setid};
      for my $val (values %$g_v) {
         my $bin = $the_histograms->bin_the_point($val);
         $the_histograms->adjust_cat_count($setid, $bin, 1);
         $min_bin = min($min_bin, $bin);
         $max_bin = max($max_bin, $bin);
      }
   }
   for my $setid (keys %sgv) { # make sure all bins have a count which is a number (not undefined).
      for my $the_bin ($min_bin .. $max_bin) { #  (0..$n_bins-1){
         $the_histograms->adjust_cat_count($setid, $the_bin, 0);
      }
   }
   $the_histograms->set_min_max_bins($min_bin, $max_bin);
   $self->{histograms} = $the_histograms;
}


# sub store_data_chunk{
#   my $self = shift;
#   my $gen_val_hashrefs = shift;
#   my @generations = keys %{$gen_val_hashrefs->[0]};

#   $self->{max_gen} = max(max(@generations), $self->{max_gen});

#   #print "ref(generations): ", ref($generations), "\n"; exit;
#   for my $i_run (0..$self->{n_runs}-1) { # 0-based
#     my $g_v = $gen_val_hashrefs->[$i_run];
#     for my $gen (@generations) {
#       my $val = $g_v->{$gen};
#       $self->store_data_point($i_run, $gen, $val);
#     }
#   }
#   if (defined $self->{histograms}) {
#     if ($self->{binnable}  and  !defined $self->{histograms}->{binning_lo_val}) { # do initial binning.
#       print "ChainData param name: ", $self->{parameter_name}, ". Doing initial binning.\n"; # sleep(1);
#       $self->bin_the_daminweighta();
#     }
#     $self->{histograms}->populate($gen_val_hashrefs);
#   } else {
#     die "**************************\n";
#   }
# }

# sub update{
#   my $self = shift;
#   my $gen_val_hashrefs = shift;
#   print stderr "In update. before store_data_chunk\n";
#   $self->store_data_chunk($gen_val_hashrefs);
#   #$self->{histograms}->populate($gen_val_hashrefs);
#   # if (defined $self->{histograms}) {
#   #    $self->{histograms}->populate($gen_val_hashrefs);
#   #  } else {
#   #    die "**************************\n";
#   #  }
#   print stderr "In update. before delete_pre_burn_in\n";

#   $self->delete_pre_burn_in();
# }


# sub get_post_burn_in_param_data_arrays{ # return array ref holding an array ref for each run, which holds param values
#   my $self = shift;
#   my @gens_big_to_small = reverse @{$self->{generations}->[0]};
#   my @runs_param_data = @{$self->{run_gen_value}}; # 
#   my @y = ();
#   my $max_gen = $gens_big_to_small[0];
#   my $end_burn_in = int($self->{burn_in_fraction} * $max_gen);
#   #print "max gen end burnin: $max_gen,  $end_burn_in \n";
#   for my $run_data (@runs_param_data) { # one element for each param
#     #print "ref(run_data): ", ref($run_data), "\n";
#     #print "run_data values: ", join(";", values %$run_data), "\n";
#     # while(my ($k, $v) = each %$run_data){
#     #   print "key, value: $k  $v  \n";
#     # }
#     my @data_array = ();
#     for my $gen (@gens_big_to_small) {
#       #  print "GENERATION: $gen  $end_burn_in.\n";
#       last if($gen <= $end_burn_in);
#       my $datum = $run_data->{$gen};
#       #  print "datum: $datum\n";
#       push @data_array, $datum;
#     }
#     #print "Data array: ", join(";", @data_array), "\n";
#     push @y, \@data_array;
#   }
#   return \@y;
# }

1;
