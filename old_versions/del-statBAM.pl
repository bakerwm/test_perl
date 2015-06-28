#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename qw(dirname basename);

my $samtools = "/home/wangming/software/samtools-1.1/samtools";
#my $ref_fai  = "/home/wangming/database/H37Rv/NC_000962.fna.fai";

my $input_dir = shift or die "Input dir that contain 02.Alignment_* directoreis, eg: ./\n";
mkdir "statBAM" unless -d "statBAM";

# Find work dir
my $work_dir = `pwd`; chomp($work_dir); $work_dir = abs_path($work_dir);

# Find 02.Alignment_* dir
my @fa = <$input_dir\/02\.Alignment\_*\/*.s.bam>;
my $num = @fa;
print "\#Find $num *.s.bam files in $input_dir\n";

die "Find no *.s.bam files in \<$input_dir\/02\.Alignment_\*\/\> \n" if(@fa < 1);

my @Run = ();
foreach my $f (@fa){
    my $f = abs_path($f);
    my $f_name = basename($f);
    my $f_dir  = dirname($f);
    my ($f_dir_sub) = $f_dir =~ /\/02\.Alignment\_(\w+)/;
    my $stat_out = $work_dir."\/statBAM\/$f_dir_sub\_$f_name\.stat";
    push @Run, "$samtools flagstat $f > $stat_out";
}
print join"\n", @Run;




