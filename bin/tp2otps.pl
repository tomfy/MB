#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw ( min max sum );
use File::Spec qw ( splitpath );
use lib '/home/tomfy/MB/lib';
use Mrbayes; # perl module encapsulating mrbayes bayesian phylogeny program.

my %newick_number = ();
my %number_newick = ();

my $tp_filename = undef;
my $format = 'long';
GetOptions(
           'tp_file=s' => \$tp_filename, # the .tp output file from running MB.pl
           'format=s' => \$format,       # long/short
          );

my ($all_lines, $last_gen) = read_tp_file($tp_filename);

my $otps_filename = $tp_filename;
$otps_filename =~ s/[.]tp//;
#print "otps filename: $otps_filename \n";
$otps_filename .= '.otps_' . (($format eq 'long')? 'long' : 'short');
open my $fhotps, ">", $otps_filename or die "Couldn't open $otps_filename for writing.\n";
#print "otps filename: $otps_filename \n";


# get number of leaves from newick expression;
my $first_line = $all_lines->[0];
my @cols = split(" ", $first_line);
my $nwck = $cols[9];
my $leaf_count = 0;
while ($nwck =~ s/\d+//) {
   $leaf_count++;
}
print "# leaf count: $leaf_count\n";

my $topo_number = 0;

for my $aline (@$all_lines) {
   my @cols = split(" ", $aline);
   my ($gen, $newick) = @cols[0,9];
   my ($minlabel, $ordered_newick, $tree_labels, $split_count) = get_ordered_newick($newick, $leaf_count);

   my $topo_string;
   if ($format eq 'long') {
      my $the_split = join( ";", sort { $a <=> $b } keys %$split_count );
      $topo_string = "$ordered_newick  $the_split"
   } else {
      if (exists $newick_number{$ordered_newick}) {
         $topo_string = $newick_number{$ordered_newick} . ':-:-';
      } else {
         $topo_number++;
         $newick_number{$ordered_newick} = $topo_number;
         my $the_split = join( ";", sort { $a <=> $b } keys %$split_count );
         $topo_string = $topo_number . ':' . $ordered_newick . ':' . $the_split;
      }
   }

   printf $fhotps ("%6d  %s %s %s %s  %3d %3d %3d   %s  %s\n", @cols[0..8], $topo_string);
}

##########################

sub get_ordered_newick{
   my $newick = shift;
   my $n_taxa = shift;
   #   $newick =~ s/;\s*$//;        # remove final semi-colon.
   $newick =~ s/^\s+//;         # remove init whitespace.
   #   $newick =~ s/:[0-9]+[.][0-9]+(e[-+][0-9]{2,3})?(,|[)])/$2/g; # remove branch lengths
   my $mask = ( 1 << $n_taxa ) - 1; # n_taxa 1's
   my $split_count = {}; # keys are split bit patterns, values: counts
   my ( $minlabel, $ordered_newick, $tree_labels, $split_bp ) =
     Mrbayes::order_newick_topo_only( $newick, $split_count, $mask, 0 );
   #  Mrbayes::order_newick( $newick, $split_count, $mask, 0 );
   $newick = $ordered_newick;
   if ( $ordered_newick =~ s/^\s*[(]1,/'(1,(/ ) { # go from form (1,(...),(...)) with trifurcation, to (1,((...),(...))) with only bifurcations
      $ordered_newick .= ")'";
   } else {
      die "newick: $ordered_newick doesn't start with '(1,' as expected.\n";
   }
   $ordered_newick =~ s/[']//g;
   return ($minlabel, $ordered_newick, $tree_labels, $split_count);
}

sub read_tp_file{
   my $tp_filename = shift;
   #   my $burn_in = shift || 0.1; # if < 1 interpret as burn-in fraction, else as max burn-in generation number.
   open my $fht, "<", "$tp_filename" or die "Couldn't open $tp_filename for reading.\n";
   my @all_lines = <$fht>;
   my $last_line = $all_lines[-1];
   my $last_gen =  ($last_line =~ /^\s*(\d+)/)? $1 : undef;
   #   my $max_burn_in_gen = ($burn_in < 1)? int($burn_in * $1) : $burn_in;
   return (\@all_lines, $last_gen);
}

