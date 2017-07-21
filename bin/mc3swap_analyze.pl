#!/usr/bin/perl -w
use strict;
use lib '/usr/local/share/perl/5.14.2/';
use List::Util qw ( min max sum );
use Math::Symbolic;
use Algorithm::CurveFit;
use Data::Dumper;
use Getopt::Long;

# analyze an .mc3swap file (produced by MB.pl), which
# has info on number of swaps attempted and accepted for each pair of T-levels;
# 

my $n_runs = undef;
my $n_temps = undef;
my $dT = undef;
my $burn_in = 0.1;
my $max_iter = 200;
GetOptions(
        #   'nruns=i' => \$n_runs, # fasta alignment file
	   'nTs=i' => \$n_temps,
           'Tstep=f' => \$dT,   # if <= 1 delta T, if > 1, Tratio
           'burnin=f' => \$burn_in, # if > 1, start with the first line with gen >= this, if < 1 fractional burn-in
          );
#print "XXX nTs, dT, burn: $n_temps $dT $burn_in\n";


my @deltats = ();
my @paccepts = ();


my @lines = <>;                 # read in all lines of mc3swap file
my $last_line = $lines[-1];
my $llcopy = $last_line;
my @xs = split(" ", $llcopy);
my $n = scalar @xs - 1;
$n_runs = $n / ($n_temps*($n_temps-1));

my @pas = ();
for (1..$n_temps) {
   push @pas, [(0) x ($n_runs+1)];
}

# print "n_runs: $n_runs nTs $n_temps n: $n \n";
my $end_burn_in_line;
if ($burn_in < 1) {
   $end_burn_in_line = $lines[ int( (scalar @lines ) * $burn_in)];
   print "# $end_burn_in_line";
} else {
   for (@lines) {
      if (/^\s*(\d+)\s/ and ($1 >= $burn_in)) {
         $end_burn_in_line = $_;
         print "# $end_burn_in_line";
         last;
      }
   }
}
print "# $last_line"; 
my @ll = split(" ", $last_line);
my @ebil = split(" ", $end_burn_in_line);
my $n_pas = ($n_temps*($n_temps-1)/2);
die "Number of mc3swap.out cols inconsistent with nruns, nTs.\n" if((scalar @ll -1) != ($n_runs * ($n_temps*($n_temps-1))));
my $lgen = shift @ll;
my $ebigen = shift @ebil;
my @pbi_cols = ();
for (0..scalar @ll -1) {
   my $diff = $ll[$_] - $ebil[$_];
   push @pbi_cols, $diff;
}
print "# ", join(" ", @pbi_cols), "\n";

my @pa_sums = ((0) x $n_pas);
my @pa_sumsqrs = ((0) x $n_pas);

while (@pbi_cols) {
   for my $i_run (1..$n_runs) {
      my $i = 0;
      for my $gap (1..$n_temps-1) {
         my $n_pairs = $n_temps - $gap; # number of temperature pairs which are $gap T-levels apart.
         for my $t_lo (0..$n_pairs-1) {
            my $t_hi = $t_lo + $gap;
            my $n_accepted = shift @pbi_cols;
            my $n_out_of = shift @pbi_cols; 
            my $p_accept = $n_accepted / $n_out_of;
            printf("%7.4f  ", $p_accept);
            $pa_sums[$i] += $p_accept;
            $pa_sumsqrs[$i] += $p_accept*$p_accept;
            if ($t_lo == 0) {
               # print "$lgen $i_run $t_lo $t_hi $p_accept \n";
               $pas[$t_hi]->[$i_run] = $p_accept;
               my $TminusTcold = ($dT <= 1.0)? $t_hi*$dT : $dT**$t_hi;
               push @deltats, $TminusTcold;
               push @paccepts, $p_accept;
            }
            $i++;
         }
         print "\n";
      }
      print "\n";
   }
}
my $i = 0;
for my $gap (1..$n_temps-1) {
my $Tratio = $dT**$gap; 
printf("# %7.4f  ", $Tratio);
   my $n_pairs = $n_temps - $gap; # number of temperature pairs which are $gap T-levels apart.
   for my $t_lo (0..$n_pairs-1) {
      my $pa_avg = $pa_sums[$i]/$n_runs;
      my $pasqr_avg = $pa_sumsqrs[$i]/$n_runs;
      my $stderr = sqrt( ($pasqr_avg - $pa_avg**2)/($n_runs-1) );
      printf("%6.4f +- %6.4f  ", $pa_avg, $stderr);
      $i++;
   }
   print "\n";
}

for my $t_hi (1..$n_temps-1) {
   my $TminusTcold = ($dT <= 1.0)? $t_hi*$dT : $dT**$t_hi;
   printf("# %3i %7.4f  ", $t_hi, $TminusTcold);
   my $sum = 0;
   for my $i_run (1..$n_runs) {
      my $the_pa = $pas[$t_hi]->[$i_run];
      $sum += $the_pa;
      printf("%7.5f ", $pas[$t_hi]->[$i_run]);
   } 
   printf("%7.5f \n", $sum/$n_runs);
}
print "# deltats: ", join(", ", @deltats), "\n";

my ($x1, $y1) = ($deltats[0], $paccepts[0]);
my ($x2, $y2) = ($deltats[-1], $paccepts[-1]);
# unshift @deltats, 1;
# unshift @paccepts, 1;
my $s = (log($y2) - log($y1))/($x2 - $x1);
my $c_guess = $y1/exp($x1*$s);
my $a_guess = $s;
print "#  initial a, c: $a_guess, $c_guess \n";
my $formula = 'c * exp(a*x)';
my $variable = 'x';
my @parameters = (['a', $a_guess, '0.000001'], ['c', $c_guess, '0.000001']);

my $square_residual = Algorithm::CurveFit->curve_fit(
                                                     formula => $formula,
                                                     params => \@parameters,
                                                     variable => $variable,
                                                     xdata => \@deltats,
                                                     ydata => \@paccepts,
                                                     maximum_iterations => $max_iter
                                                    );
my $a_best = $parameters[0]->[1];
my $c_best = $parameters[1]->[1];
print "#  best a, c: $a_best, $c_best \n";
print "#  rms residual: ", ($square_residual/($n_runs*($n_temps-1)))**0.5, "\n";

my ($dt_best, $pa_best, $max_mstj) = (-1, -1, -1);
$dt_best = -2/$a_best;
$pa_best = $c_best * exp($a_best * $dt_best);
$max_mstj = $pa_best * $dt_best**2;
printf("#  Opt. Tgap, Pa, msTjump: %8.5f %8.5f %8.5f \n", $dt_best, $pa_best, $max_mstj);
my $pr_swap = 0.234;
my $Trat_prswap = log($pr_swap/$c_best)/$a_best;
print "# Use Tratio of:  $Trat_prswap  for swap prob of:  $pr_swap \n";
