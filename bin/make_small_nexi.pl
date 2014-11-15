#!/usr/bin/perl -w
use strict;

my $whole_nexus_file = shift;

my $chunk_size = shift || 100000;
my $n_taxa = shift || 5;

open my $fh, "<", "$whole_nexus_file" or die "Couldn't open $whole_nexus_file for reading.\n";

my $top_string = '';
my @alignment_lines;

while (<$fh>) {
  $top_string .= $_;
  last if(/^\s*matrix/);
}
$top_string =~ s/nchar=(\d+)/nchar=$chunk_size/;
while (<$fh>) {
  last if(/^\s*;/);
  push @alignment_lines, $_;
}

my $middle_string = '';
my $bottom_string =
  "; \n" .
  "endblock; \n" .
  "begin mrbayes; \n" .
  "Outgroup potato; \n" .
  "mcmc Mcmcdiagn=yes ngen=200000 samplefreq=100 printfreq=400 diagnfreq=800 relburnin=yes burninfrac=0.2 nruns=4 nchains=2 temp=0.5; \n" .
  "sumt burnin=250 calctreeprobs=yes showtreeprobs=yes; \n" .
  "lset applyto=(all) nst=6 rates=invgamma; \n" .
  "prset applyto=(all) ratepr=variable topologypr=uniform ; \n" .
  "end; \n";


my $chunk_count = 1;
while (1) {
  my $character_count = 0;
  while ($character_count < $chunk_size and scalar @alignment_lines >= 5) {
    for (1..($n_taxa+2)) {
      my $next_line = shift @alignment_lines;
      $middle_string .= $next_line;
    }
    $character_count += 50;
  }
  my $chunk_nexus_filename = "chunk$chunk_count.nexus";
  open my $fh_out, ">", "$chunk_nexus_filename";
  print $fh_out "$top_string \n" . "$middle_string \n" . "$bottom_string \n";
close $fh_out;
  $middle_string = '';
  $chunk_count++;
# last if($chunk_count > 2);
  last if(scalar @alignment_lines * 50  <  $chunk_size*7);
}



