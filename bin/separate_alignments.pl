#!/usr/bin/perl -w
use strict;
use Cwd;

my $count = 1;
while(1){
	my $nexus = read_one_nexus();
	print $nexus;
	last if(!$nexus);
	my $new_dir_name = $count;
	print "new dir name: $new_dir_name\n";
	if(mkdir $new_dir_name){
		chdir $new_dir_name or die "Couldn't change dir to $new_dir_name.\n";
		print "working dir: ", getcwd(), "\n";
	my $outfilename = 'sim_alignment.nexus';
		open my $fh, ">", "$outfilename" or die "Couldn't open file $outfilename for writing.\n";
		print $fh $nexus;
		close $fh;
		chdir '../';
	}else{
		die "couldn't make dir $new_dir_name.\n"
	}
$count++;
}


sub read_one_nexus{
	my $nexus = '';
	while(<>){
		$nexus .= $_;
# print $nexus, "\n";
		last if(/END;/);
	}
	print $nexus, "\n";
	return ($nexus =~ /[#]NEXUS/m)? $nexus : '';
}
