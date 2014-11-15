#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my %number_id_hash = ();
my $filename_base = shift;
my $requested_gen = undef;
my $run = 1;
my $translate = 1; # whether to translate number identifiers into id strings (e.g. AT1g034567.1 )

GetOptions(
    'gen=i'           => \$requested_gen,
    'run=i'     => \$run, 
    'translate!' => \$translate, # exclamation point means can use  -translate , -notranslate
);
my $tree_filename = "$filename_base.run$run.t";
my $param_filename =  "$filename_base.run$run.p";
open my $fh, "<", "$tree_filename" or die "couldnt open $tree_filename for reading.\n";

my $the_newick = '';
while (<$fh>) {
  if (/^\s*begin trees;/) {
    my $next_line = <$fh>;
      if ($next_line =~ /^\s*translate/) {
	while (<$fh>) {
	  if (/^\s*(\d+)\s+(\S+)[,;]\s*$/) {
#	   print "XXX: $1  $2  \n";
	    $number_id_hash{$1} = $2;
	  } elsif (/^\s*tree\s+gen[.](\d+)\s*[=]\s*\S+\s+([(].*[)]);/) {
	    my ($gen, $newick) = ($1, $2);
	    if ($gen == $requested_gen) {
	      if($translate){
	      while (my ($n, $id) = each %number_id_hash) {
	#	print "$n  $id \n";
		$newick =~ s/([(,])$n[:]/$1$id:/;
	#	print "$newick \n";
	      }
	    }
	      $the_newick = $newick;
	     # print $newick, "\n";
	    }
	  }
	}
      }
  }
}

close $fh;

open $fh, "<", "$param_filename" or die "couldnt open $param_filename for reading.\n";
my ($gen, $lnl);

while(<$fh>){
next if(/^\s*\D/); # skip if first non-whitespace is not digit
  my @cols = split(" ", $_);
($gen, $lnl) = @cols[0..1];
last if($gen == $requested_gen);
}

print "$requested_gen \n";
print "$the_newick;\n";
printf("%10.1f \n", $lnl);

