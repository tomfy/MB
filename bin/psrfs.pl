#!/usr/bin/perl -w
use strict;
use lib '/home/tomfy/myperlmodules/';
use Misc qw ( mean_and_var  psrf );
use Getopt::Long;
# use List::Util qw ( min max sum );



my $dir = undef;
my $burn_in_fraction = 0.25;
my $initial_group_number = 1;
GetOptions(
    'dir=s'           => \$dir,
    'init_group=s' => \$initial_group_number,
    'burn_in_fraction=f' => \$burn_in_fraction,
);
my @subdirs;
if(defined $dir){
@subdirs = ($dir);
}else{
@subdirs = split(" ", `ls`);
}

for my $adir (@subdirs){
next if(! -d $adir);
$adir =~ s/\s+$//;
print "# directory: [$adir]\n";
if($adir =~ /\/$/){
}else{
$adir .= '/';
}

# print "dir $dir \n";
my %grun_filename = ();
my @files = `ls  $adir*.run*.p`;
# print join("; ", @files), "\n";
my ($i_group, $i_run) = (-1, -1);
for (@files) {
  s/\s+$//;
	# print "file: $_\n";
  if ( /(\d+)[.].*run(\d+)[.]p/) {
    ($i_group, $i_run) = ($1, $2);
	# print "$i_group  $i_run \n";
  }
  $grun_filename{"$i_group:$i_run"} = $_;
}

$i_group = $initial_group_number-1;
while (1) {
  $i_group++;
  $i_run = 1;
	last if ($i_group > 10000);
  my $filename;
  if (exists $grun_filename{"$i_group:$i_run"}) {
    $filename = $grun_filename{"$i_group:$i_run"};
  } else {
    next;
  }
  #  last if( ! -f $filename);
  #  print "$i_group  $i_run  FILENAME [$filename] \n";
  #exit;
  my @lnl_means = (); my @lnl_vars = ();
  my @tl_means = (); my @tl_vars = ();
# print "filename: $filename \n";
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

    if (exists  $grun_filename{"$i_group:$i_run"}) {
      $filename = $grun_filename{"$i_group:$i_run"};
    } else {
      my ($lnl_psrf, $nlnl, $Blnl, $Wlnl) = psrf(scalar @trunc_lnls, \@lnl_means, \@lnl_vars);
      my ($tl_psrf, $ntl, $Btl, $Wtl) = psrf(scalar @trunc_tls, \@tl_means, \@tl_vars);
      print "$i_group  $nlnl $lnl_psrf $Blnl $Wlnl    $ntl $tl_psrf $Btl $Wtl\n";
      last;
    }
  }				# loop over runs
}				# loop over groups

}
