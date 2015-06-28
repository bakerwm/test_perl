#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Std;

################################
# filt the sRNA candidates
# 1. keep RNAz hits
# 2. TPM > cutoff
#
################################

sub usage {
    die(qq/
Usage: filt_tags.pl [options] infile.txt

Options: -e   : The column of TPM value
         -z   : The column of RNAz value
         -c   : The cutoff of TPM value for each seq [10]

Example: 
1. Filter the file by TPM > 10:

filt_tags.pl -e 14 -z 21 -c 10 Lib01.report.txt
\n/);
}


usage() if(@ARGV == 0);
my %opts = (e => 14, z => 21, c => 10);
getopts("e:z:c:", \%opts);

my $infile = shift(@ARGV);
my $notinclude = $infile . '.del';
open my $fh_in,  "< $infile" or die "Cannot open $infile, $!\n";
open my $fh_del, "> $notinclude" or die "Cannot open $notinclude, $!\n";
my $e = $opts{e} - 1;
my $z = $opts{z} - 1;
my $count_col = $e - 1;

while(<$fh_in>) {
    chomp;
    my $line = $_;
    my @ps = split(/\s+/, $_);
    my @infs = split(/\s+/, $line, 13);
    pop(@infs);
    next if(/(^\s*$)|(^\#)/);
    if($ps[$z] eq '-' ) {
        if($ps[$e] > $opts{c}) {
            print join("\t", @infs, $ps[$count_col], $ps[$e], $ps[$z]) . "\n";
        }else {
            print $fh_del $_ . "\n";
        }
    }else {
        print join("\t", @infs, $ps[$count_col], $ps[$e], $ps[$z]) . "\n";
    }
}
close $fh_in;
close $fh_del;

