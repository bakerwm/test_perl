#!/usr/bin/perl  -w
#use Time::Local;
use strict;
use POSIX qw(strftime); # For print present time.

my $usage   = 'perl  run_statPosCov.pl  match_genome_dir  '."\n";

my $dir     = shift or die $usage;

my $in_file = 'match_genome';
my @match   = glob"$dir$in_file*txt";

print "Calculate reads/5'End/3'End coverage\:\n";
if(!@match){
    print "Warn:
1. Copy match_genome files in  $dir to current work directory,
2. Change file names to \"match_genome.200PE.txt\".\n";
    exit (1);
}

my $stat_pos = '/home/wangming/work/bin/StatCovrage.pl';

foreach(@match){
    system"perl  $stat_pos  -i  $_ ";

    print strftime("%Y-%m-%d %H:%M:%S", localtime(time));
    print "\t$_\tDone!\n";
}

print strftime("%Y-%m-%d %H:%M:%S", localtime(time));
print "\tFinish!\n";
