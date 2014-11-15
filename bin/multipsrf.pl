#!/usr/bin/perl
use strict;

my $prefix = shift || 'chunk';
my $burninfrac = shift || 0.25;
my $every = shift || 1;

for(1..844){
my $filename = "$prefix$_.nexus";
last if(! -f "$filename.run1.p");
print `psrf.pl $filename $burninfrac  $every`;
}


