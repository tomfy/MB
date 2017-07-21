#!/usr/bin/perl -w
use strict;
use Cwd;
use Getopt::Long;
#my $n_runs = undef;
my $n_temps = undef;
my $dT = undef;
my $burn_in = 0.1;
my $max_iter = 200;
my $align_filename = undef;
my %seqlength_pastr = ();
GetOptions(
           #  'nruns=i' => \$n_runs, # fasta alignment file
	   'nTs=i' => \$n_temps,
           'Tstep=f' => \$dT,   # if <= 1 delta T, if > 1, Tratio
           'burnin=f' => \$burn_in, # if > 1, start with the first line with gen >= this, if < 1 fractional burn-in
           'align_filename=s' => \$align_filename,
          );

print "# $n_temps $dT  $burn_in  $align_filename \n";

my $dirs = get_directories();

my $mc3swanout_filename = 'zzzz';

print "#  ", join("; ", @$dirs), "\n";

for my $adir (@$dirs) {
   chdir($adir);
   my $cwd =  getcwd();
   #  print "cwd: $cwd \n";
   my $seqlength = undef;
   if ($cwd =~ /M(\d+)\s*$/) {
      $seqlength = $1;
   }
   my $rep_dirs = get_directories();
   my @pa_sums = ();
   my @pa_sumsqrs = ();
   my $count = 0;
   my $n_pas = 0;
   for my $rep_dir (@$rep_dirs) {
      chdir($rep_dir);
      if (-f $align_filename) {
         my $command_string = "mc3swap_analyze.pl -nTs $n_temps -Tst $dT -burn $burn_in  < $align_filename > $mc3swanout_filename";
     #    print $command_string, "\n";
         system "$command_string"; # "mc3swap_analyze.pl -nTs $n_temps -Tst $dT -burn $burn_in  < $align_filename > $mc3swanout_filename";
#print "YYY\n";
      }
      if (-f $mc3swanout_filename) {
         my ($count_one_file, $pa_sums_one_file) = read_zzz($mc3swanout_filename);
      #   print "$count_one_file  ", join(", ", @$pa_sums_one_file), "\n";
         $count += $count_one_file;
         while (my ($j, $pa) = each @$pa_sums_one_file) {
            $pa_sums[$j] += $pa;
            $pa_sumsqrs[$j] += $pa*$pa;
         }
      }
      chdir('../');
   }
   my $str = sprintf("%4i   ", $count);
   while (my($j, $pa_sum) = each @pa_sums) {
      my $pa_avg = $pa_sum/$count;
      my $pasqr_avg = $pa_sumsqrs[$j]/$count;
      my $stderr = sqrt( ($pasqr_avg - $pa_avg**2)/($count*($count-1)) );
      $str .= sprintf("%8.5f +- %8.5f  ", $pa_avg, $stderr);
   }
   $seqlength_pastr{$seqlength} = $str;
   chdir('../');
}
my @sseqlengths = sort {$a <=> $b} keys %seqlength_pastr;

for (@sseqlengths) {
   printf("%5i %s \n", $_, $seqlength_pastr{$_});
}

sub get_directories{
   my @lsoutlines = `ls -l`;
   my @dirs = ();
   for (@lsoutlines) {
      if (/^d.*\s(\S+)\s*$/) {
         push @dirs, $1;
      }
   }
   # print join("\n", @dirs), "\n";
   return \@dirs;
}

sub read_zzz{
   my $filename = shift;
   open my $fh, "<", "$filename" or die "couldn't open $filename for reading.\n";
   my @pas = ();
   my ($count, $i) = (0, 0);
   while (<$fh>) {
      next if(/^\s*#/);
      if (/^\s*$/) {
         $i = 0;
         $count++;
      } else {
         my @xs = split(" ", $_);
         for my $x (@xs) {
            $pas[$i] += $x;
            $i++;
         }

      }
   }
   #   print "$count  \n";
   #   print join(" ", @pas), "\n";
   return ($count, \@pas);
}
