#!/usr/bin/perl -w
use strict;

my @nexus_files = `ls *.nexus`;

if(scalar @nexus_files == 0){
my @run1t_files = `ls *.run1.t`;
for(@run1t_files){
	s/[.]run1[.]t\s*$//; 
	push @nexus_files, $_;
}
}

print join("\n", @nexus_files), "\n";
# exit;

for my $nexus_file (@nexus_files) {
  chomp $nexus_file;

  my $psrf_output = `psrf.pl $nexus_file 2> /dev/null `;
next if($psrf_output =~ /(Only one |No )run[*][.][tp]/);
  my @psrf_cols = split(" ", $psrf_output);


  my $topolyzer_output = `topolyzer.pl $nexus_file 2> /dev/null `;
print $topolyzer_output, "\n";
exit;
next if($topolyzer_output =~ /(Only one |No )run[*][.][tp]/);
  my @topo_lines = split("\n", $topolyzer_output);
  my $topo_max_tvd = -1;
my $topo_clstr_tvd = -1; # divide runs into 2 groups, avg tvd between 2 groups, maximized over possible ways of dividing runs into 2 groups.
  for (@topo_lines) {
    if (/^\s*#\s*topology[.]\s*L1:\s*\(.*\s+([\d.]+)\s*\)\s*([\d.]+)/) {
      $topo_max_tvd = $1;
      $topo_clstr_tvd = $2;
      last;
    }
  }

  print $nexus_file, "   ", $psrf_cols[1], "  ", $psrf_cols[5], "  ", $topo_max_tvd, "  ", $topo_clstr_tvd, "\n";


}
