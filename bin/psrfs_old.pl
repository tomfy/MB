#!/usr/bin/perl -w
use strict;
use List::Util qw ( min max sum );

my @subdirs;
my $dir = shift;
if(defined $dir){

@subdirs = ($dir);
}else{
@subdirs = split(" ", `ls`);
}

for my $adir (@subdirs){

$adir =~ s/\s+$//;
print "# directory: [$adir]\n";

if($adir =~ /\/$/){
}else{ 
$adir .= '/'; 
}

# print "dir $dir \n";
my %grun_filename = ();
my @files = `ls  $adir*group*.*.run*.p`;
my ($i_group, $i_run) = (-1, -1);
for (@files) {
  s/\s+$//;
  if ( /group(\d+)[.].*run(\d+)[.]p/) {
    ($i_group, $i_run) = ($1, $2);
  }

  $grun_filename{"$i_group:$i_run"} = $_;
}

$i_group = 0;

#print join("\n", @files), "\n";
#exit;
while (1) {
  $i_group++;
  $i_run = 1;
  my $filename;
  if (exists $grun_filename{"$i_group:$i_run"}) {
    $filename = $grun_filename{"$i_group:$i_run"};
  } else {
    last;
  }
  #  last if( ! -f $filename);
  #  print "$i_group  $i_run  FILENAME [$filename] \n";
  #exit;
  my @lnl_means = (); my @lnl_vars = ();
  my @tl_means = (); my @tl_vars = ();

  while ( -f $filename ) {	# loop over runs
    my @lnls = ();
    my @tls = ();
    open my $fh, "<", "$filename";
    <$fh>; <$fh>;		#
    while (<$fh>) {
      #   print "XXXXX $_ \n";
      my @cols = split(" ", $_);
      my ($gen, $lnl, $tl) = @cols[0..2]; # ($cols[0], $cols[1], $cols[2]); #
      push @lnls, $lnl;
      push @tls, $tl;
      #   print "$lnl, $tl \n";
    } # loop over lines in file
    #exit;
    my ($mean_lnl, $var_lnl) = mean_and_var(\@lnls);
    my ($mean_tl, $var_tl) = mean_and_var(\@tls);
    #print "XXX $mean_lnl, $var_lnl, $mean_tl, $var_tl \n";
    push @lnl_means, $mean_lnl;
    push @lnl_vars , $var_lnl;
    push @tl_means, $mean_tl;
    push @tl_vars , $var_tl;

    $i_run++;
    
    if (exists  $grun_filename{"$i_group:$i_run"}) {
      $filename = $grun_filename{"$i_group:$i_run"};
    } else {
      my $lnl_psrf = psrf(scalar @lnls, \@lnl_means, \@lnl_vars);
      my $tl_psrf = psrf(scalar @tls, \@tl_means, \@tl_vars);
      print "$i_group  $lnl_psrf  $tl_psrf \n";
      last;
    }
  }				# loop over runs
}				# loop over groups

}
sub mean_and_var{
  my $xs = shift;		# array ref
  my $burn_in_fraction = shift;
  $burn_in_fraction = 0.25 if( !defined $burn_in_fraction);
  my ($x_sum, $xsq_sum) = (0, 0);
  my $n = scalar @$xs;
  if ($n > 1) {
    my $n_start = int($burn_in_fraction * $n);
    #print "$n_start $n \n";
    #print join("\nzzz ", @$xs), "\n";
    #exit;
   
 
    my $n_in_sum = $n - $n_start;
    my @xarray = @$xs;
    my $mean = sum(@xarray[$n_start .. $n-1])/$n_in_sum;
    my $var = 0;
    for (my $i = $n_start; $i < $n; $i++) {
   #   my $x = $xs->[$i];
      #   $x_sum += $x;
   #    $xsq_sum += $x**2; # $x*$x;
     $var += ($xs->[$i] - $mean)**2;
    }
    $var /= ($n_in_sum - 1);
    #  my $mean = $x_sum/$n_in_sum;
 #   my $var1 = ($xsq_sum/$n_in_sum - $mean*$mean)*$n_in_sum/($n_in_sum - 1);;
#    print "n, mean, var: $n   $n_in_sum $mean  $var \n";
    return ($mean, $var);
  } elsif ($n == 1) {
    return ($xs->[0], undef);
  } else {
    return (undef, undef);
  }
}

sub psrf{
  my $n = shift;		# values in each mean
  my $means = shift;
  my $vars = shift;
  my $m = scalar @$means;
  my $mean_of_means = sum(@$means)/$m;
  my $W = sum(@$vars)/$m;
#print "X  $mean_of_means, $mean_of_vars \n";
  my $B = 0;
  for my $chain_mean (@$means) {
    $B += ($chain_mean - $mean_of_means)**2;
  }
  $B /= ($m-1);

#print "XX W, B: $W, $B \n";
  return sqrt(1 - 1/$n + $B/$W);
}
