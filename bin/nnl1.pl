#!/usr/bin/perl -w
use strict;

my $prev = [];
my $this = [];
while(<>){
	next if(/^\s*#/); # skip comment lines
		my @cols = split(" ", $_);
	my $chunk_number = shift @cols;
	my $total = pop @cols;
	$this = \@cols;
	my $L1 = 0;
#	print "prev: ", join(", ", @$prev), " this: ", join(", ", @$this ), "\n";
	if(scalar @$prev == scalar @$this){
		for my $i (0 .. scalar @$this -1){
			$L1 += abs($prev->[$i] - $this->[$i]);
		}
		$L1 /= 2*$total;
		print "$chunk_number  $L1 \n";
	}		
	$prev = $this;
}
