#!/usr/bin/perl -w
use strict;

# parse output from GeoMeTree
# which gives the 'geodesic' path 
# between two trees in terms of the set
# of splits (and associated weights aka branch lengths)
# along the path.
# first we need to get all the splits (specified in some consistent way)
# and their names.


my %id_dummytax = ();
my $name_split = {};
while(<>){

# get dummy taxa:
	if(/dummy taxa/){
		while(1){
			my $line = <>;
			# print "$line \n";
			if($line =~ /^\s*(\S+)\s+(\S+)/){
				$id_dummytax{$1} = $2;
			}else{
				last;
			}
		}
		while(my ($i, $t) = each %id_dummytax){
		#	print "dummy id, taxon: $i  $t \n"; 
		}	
	}
	if(/^\s* Splits \s (only \s in | common \s to \s both) /x){
		get_splits($name_split);
	}

}

for(keys %$name_split){
#	print "$_  ", $name_split->{$_}, "\n";
	print "$_  ", format_split($name_split->{$_}, \%id_dummytax), "\n";
}

sub get_splits{

	my $name_split = shift;

	while(<>){

		if(/^ \s* (\d+ \/ \d+) \s+ (\S+) \s+ ([0-9.]+) \s+ ([0-9.]+) \s* $/x ){
			$name_split->{$1} = $2;
		}else{
			last;
		}
	}
	return $name_split;
}

sub format_split{
	my $split_string = shift;
	my $dummy_id_tax = shift;

	while ( my ($id, $name) = each %$dummy_id_tax){
	#	print "XXX: $id  $name    ";
		$split_string =~ s/$id/$name/; # replace dummy ids with numbers
		#	print "A: $split_string \n";
	}

	my @ids = split("[*]", $split_string);
#	@ids = sort {$a <=> $b} @ids;
	$split_string = join(",", @ids);
	return $split_string;
}	


