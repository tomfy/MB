#!/usr/bin/perl -w
use strict;
use List::Util qw( sum min max);
use lib '/home/tomfy/myperlmodules';
use Misc qw ( mean_and_var  psrf );

my $filenamebase = shift; # e.g. group1.nexus
my $burn_in_fraction = shift || 0.25;
my $every = shift || 1;
# print "dir $dir \n";
my %grun_filename = ();
my @files = `ls  $filenamebase.run*.p`;
if(scalar @files == 0){
  print "$filenamebase  No run*.p files.\n";
  exit;
}elsif(scalar @files == 1){
  print "$filenamebase Only one run*.p file. No psrf info.\n";
exit;
}
my $i_run = -1;
for (@files) {
  s/\s+$//;
  if ( /[.]run(\d+)[.]p/) {
    $i_run = $1;
  }

  $grun_filename{"$i_run"} = $_;
}

#while (1) {
$i_run = 1;
my $filename;
if (exists $grun_filename{"$i_run"}) {
  $filename = $grun_filename{"$i_run"};
} else {
  last;
}
my @lnl_means = (); my @lnl_vars = ();
my @tl_means = (); my @tl_vars = ();


while ( -f $filename ) {	# loop over runs
  my @lnls = ();
  my @tls = ();
  open my $fh, "<", "$filename";
  $_ = <$fh>; 
# print "$filename  $_";
$_ = <$fh>;			#
#	print "$filename $_ \n";
my $line_count = 0;  
while (<$fh>) {
	$line_count++;
	next if($line_count % $every  != 0);
    #   print "XXXXX $_ \n";
    my @cols = split(" ", $_);
    my ($gen, $lnl, $tl) = @cols[0..2]; # ($cols[0], $cols[1], $cols[2]); #
    push @lnls, $lnl;
    push @tls, $tl;
    #   print "$lnl, $tl \n";
  } # loop over lines in file
    #exit;

	my $n_start = int($burn_in_fraction* scalar @lnls);
my @trunc_lnls = @lnls[$n_start..scalar @lnls -1];
my @trunc_tls = @tls[$n_start..scalar @tls -1];
  my ($mean_lnl, $var_lnl) = mean_and_var(\@trunc_lnls);
  my ($mean_tl, $var_tl) = mean_and_var(\@trunc_tls);
  #print "XXX $mean_lnl, $var_lnl, $mean_tl, $var_tl \n";
  push @lnl_means, $mean_lnl;
  push @lnl_vars , $var_lnl;
  push @tl_means, $mean_tl;
  push @tl_vars , $var_tl;

  $i_run++;

  if (exists  $grun_filename{"$i_run"}) {
    $filename = $grun_filename{"$i_run"};
  } else {
    my ($lnl_psrf, $nlnl, $Blnl, $Wlnl) = psrf(scalar @trunc_lnls, \@lnl_means, \@lnl_vars);
    my ($tl_psrf, $ntl, $Btl, $Wtl) = psrf(scalar @trunc_tls, \@tl_means, \@tl_vars);
    print "$filenamebase  $lnl_psrf  $nlnl $Blnl $Wlnl    $tl_psrf $ntl $Btl $Wtl \n";
    last;
  }
}
