#!/usr/bin/perl -w
use strict;
use lib '/usr/share/perl/5.14.2/';
use Graphics::GnuplotIF qw(GnuplotIF);
my $input_filename = shift; # 

my $persist = 0;
my $plot1 = Graphics::GnuplotIF->new( persist => $persist, style => 'lines lw 2');

#$plot1->gnuplot_cmd(' set style data histograms ');
#$plot1->gnuplot_cmd(' set style histogram rowstacked ');
#$plot1->gnuplot_cmd(' set style fill solid 1.0 border 1 ');

#$plot1->gnuplot_cmd(' set origin 0.0, 0.75 ');
# $plot1->gnuplot_cmd(' plot [0:6] sin(x) ');
$plot1->gnuplot_cmd(" plot [][0:6400] \
     'topology_pd_vs_chunks_1-171' using 1:3 t'' "); #, \
#      '' using 1:2 t'', \
#      '' using 1:3 t''
# ");

 $plot1->gnuplot_pause(0);


# my @lines = <$fh>;
# while(1){
# my $a_line = shift @lines;
# if($a_line =~ /^\s*#/){}
# else{ unshift @lines, $a_line; last}
# }

# set style data  histograms
# set style histogram rowstacked
# set style fill solid 1.0 border 1
