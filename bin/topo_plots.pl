#!/usr/bin/perl -w
use strict;

my $input_filename = shift;


my $infilename = '"filename=' . "'" . "$input_filename" . "'" . '"';
my $outfilename = '"outfilename=' . "'" . "$input_filename.png" . "'" . '"';
my $command = 'gnuplot -e ' . $infilename . 
  # "  " . $outfilename . 
  ' ~/MB/bin/topo_plots.gnuplot';
system "$command";
