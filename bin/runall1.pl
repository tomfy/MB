#!/usr/bin/perl -w
use strict;

my $first = shift || 1;
my $last = shift || 3;
my $prefix = shift || 'chunk'; 

print STDERR "$prefix $first $last";
my $chunks_done = 0;
for my $number ($first .. $last){
my $nexus_name = "$prefix$number.nexus";
if(! -f $nexus_name){
warn "$nexus_name is not a regular file.\n";
next;
}	
my $stdout_string = `mb $nexus_name `;
	my $stdout_filename = "$prefix$number.stdout";
	print STDERR "stdout filename: $stdout_filename \n";

	open my $fh, ">", "$stdout_filename"; # or die "couldn't open $stdout_filename for writing \n";
	print $fh $stdout_string;
	close $fh;
	print STDERR "Done with chunk $number.\n";
	$chunks_done++;
	print STDERR "$chunks_done chunks done.\n";
	last if($chunks_done > 10000);
}



exit;
