package Cluster;
use strict;
use List::Util qw ( min max sum );

sub new {
  my $class             = shift;
  my $labels = shift; # arrayref holding the labels of the items being clustered
  my $lp_d = shift; # a hashref. keys index pairs (e.g. '0,2' ), values distances
  my $self = bless {}, $class;
$self->{cluster_names} = {0 => '0', 1 => '1', 2 => '2', 3 => '3',
		     4 => '4', 5 => '5', 6 => '6', 7 => '7',
		     8 => '8', 9 => '9', 10 => 'A', 11 => 'B',
		     12 => 'C', 13 => 'D', 14 => 'E', 15 => 'F'};

  my %label_cl = ($labels->[0] => '0');
  for (my $i = 1; $i < scalar @$labels; $i++) {
    $label_cl{$labels->[$i]} = '0';
  }
 
  $self->{label_cluster} = \%label_cl; 
  $self->{labels} = $labels;
$self->{labelpair_distance} = $lp_d;
$self->{label_cluster_best} = {};
  return $self;
}

sub _next{
  my $self = shift;
  my $n_clusters = shift;
  my $labels = $self->{labels};
  my $label_cluster = $self->{label_cluster};
  my $index = 1;
  my $max_allowed_digit = 1;
  # my $max_digit = $max_lower_order_digit + 1;
  while (1) {
    my $the_label = $labels->[$index];
    if ( $label_cluster->{$the_label} < $max_allowed_digit) {
      $label_cluster->{$the_label} += 1;
      last;
    # } elsif ($label_cluster->{$the_label} == $max_lower_order_digit) {
    #   $label_cluster->{$the_label} += 1;
    #   $index++;
    #   $max_lower_order_digit++;
    #   last if($index == scalar @$labels);
    } elsif ($label_cluster->{$the_label} == $max_allowed_digit) {

      $label_cluster->{$the_label} = 0;
      $index++;
 $max_allowed_digit = min($n_clusters-1, $max_allowed_digit+1);
      last if($index == scalar @$labels);
    }
  }
      return $self->partition_as_string;
}

sub next{
my $self = shift;
my $n_clusters = shift || 2;

for(0..50){
  my $partition_string = $self->_next($n_clusters);
#print "  $partition_string \n";
return $partition_string if(check_partition_string($partition_string));
}
}

sub check_partition_string{
my $string = shift;
my @p = reverse( split("", $string));
my $big_dig_so_far = 0;
for(@p){
 # print "     $big_dig_so_far   $_ \n";
  return 0 if($_ > $big_dig_so_far+1);
    $big_dig_so_far = max($big_dig_so_far, $_);
}
return 1;
}

sub partition_goodness{
  my $self = shift;

  my $labels = $self->{labels};
  my $label_cluster = $self->{label_cluster};
  my ($count_inter, $sum_d_inter, $count_intra, $sum_d_intra) = (0, 0, 0, 0);
my @inter_dists; my @intra_dists;
  for (my $i = 0; $i < scalar @$labels; $i++) {
    my $l1 = $labels->[$i];
    for (my $j = $i+1; $j < scalar @$labels; $j++) {
      my $l2 = $labels->[$j];
      my $lp = "$l1,$l2";
      my $distance = $self->{labelpair_distance}->{$lp};
      if ($label_cluster->{$l1} eq $label_cluster->{$l2}) { # intracluster
	$sum_d_intra += $distance; $count_intra++;
	push @intra_dists, $distance;
      } else { # intercluster
	$sum_d_inter += $distance; $count_inter++;
	push @inter_dists, $distance;
  }
    }
  }
  my ($avg_intra, $avg_inter) = (0, 0);
# in count is zero, sum will be zero, set avg to 0
  if ($count_intra > 0){
    $avg_intra = $sum_d_intra/$count_intra;
  }
  if($count_inter > 0){
    $avg_inter = $sum_d_inter/$count_inter;
  }
  return ($avg_inter, $avg_intra, \@inter_dists, \@intra_dists);
}

sub partition_as_string{
  my $self = shift;
  my $string = '';
  for (@{$self->{labels}}) {
    $string = $self->{cluster_names}->{$self->{label_cluster}->{$_}} . $string;
  }
  return $string;
}

sub best_n_partition{
  my $self = shift;
  my $n_clusters = shift || 2;
my $criterion = shift || 'max_inter'; #
  my $init_bip_string =  $self->partition_as_string();
  my ($best_avg_d_inter, $best_avg_d_intra) = $self->partition_goodness();
  my $bestness = $best_avg_d_inter;
  if($criterion eq 'max_inter_min_intra'){ $bestness -= $best_avg_d_intra; }
  my $best_bip_string = $init_bip_string;
  #print "$init_bip_string   $bestness  $best_avg_d_inter  $best_avg_d_intra \n";
my %partition_score = ();
my $bip_string = $init_bip_string;
my ($inter_dists, $intra_dists, $best_inter_dists, $best_intra_dists) = ([], [], [], []);
  while (1) {
 
    my ($avg_d_inter, $avg_d_intra, $inter_dists, $intra_dists)  = $self->partition_goodness();
    my $goodness = $avg_d_inter;
    $goodness -= $avg_d_intra if($criterion eq 'max_inter_min_intra');;
    $partition_score{$bip_string} = $goodness;
#   print "AAAXXX: $bip_string   $goodness  $avg_d_inter  $avg_d_intra \n";
    if ($goodness > $bestness) {
      $bestness = $goodness;
      $best_bip_string = $bip_string;
      $best_avg_d_inter = $avg_d_inter;
      $best_avg_d_intra = $avg_d_intra;
      $best_inter_dists = $inter_dists;
  $best_intra_dists = $intra_dists;
      for (keys %{$self->{label_cluster}}) {
	$self->{label_cluster_best}->{$_} = $self->{label_cluster}->{$_};
      }
    #  print "# best so far: $best_bip_string   $bestness  $best_avg_d_inter  $best_avg_d_intra \n" ;
    }
   $bip_string =  $self->next($n_clusters);
#    print STDERR "AXX: new bip: $bip_string $init_bip_string \n";
    last if($bip_string eq $init_bip_string);
  }
  # my @sorted_partitions = sort { $partition_score{$b} <=> $partition_score{$a} } keys %partition_score;
  # for(@sorted_partitions[0..5]){
    
  #   print "#  $_  ", $partition_score{$_}, "\n";
  # }
# print "intra, inter dists: ", join(", ", sort {$a <=> $b} @$best_intra_dists), "    ",  join(", ", sort {$a <=> $b} @$best_inter_dists), "\n";
  return ($bestness, "# Best partition: $best_bip_string   $bestness  $best_avg_d_inter  $best_avg_d_intra \n");
}

1;
