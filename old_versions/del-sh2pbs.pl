#!/usr/bin/perl -w
use warnings;
use strict;
use Cwd 'abs_path';
use POSIX qw(strftime);

my $sh_dir = shift or die "Input the dir of *.sh files \n";
my @sh_files = <$sh_dir\/*.sh>;
die "No *.sh files in $sh_dir .\n" if(@sh_files < 1);

my $pbs_file = shift or die "Input the name of output pbs commands: eg: out.job";

my $work_dir = `pwd`; chomp($work_dir);
my $date = strftime "%Y%m%d", localtime;
my $pbs_name = "$work_dir\/$date\_$pbs_file";

my $pbs_header = "\#PBS -S /bin/bash \n#PBS -q Debug";
my $pbs_err    = "\#PBS -e $pbs_name\.err";
my $pbs_out    = "\#PBS -o $pbs_name\.out";
my $pbs_line   = "cd \$PBS_O_WORKDIR";

my @pbs_sh;
foreach my $sh (@sh_files){
    my $sh = abs_path($sh);
    push @pbs_sh, "sh $sh > $sh\.log ";
}
my $pbs_cmd = join"\n", ($pbs_header, $pbs_err, $pbs_out, $pbs_line, @pbs_sh);

open OUT, ">$pbs_name" or die;
print OUT $pbs_cmd,"\n";
close OUT;

print "Write the following commdands to: $pbs_name\n\n";

print $pbs_cmd,"\n";

