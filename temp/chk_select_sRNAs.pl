#!/usr/bin/env perl

#########################################
# Filt the sRNAs by :
# 1. length (-l)
# 2. TPM (-e)
# 3. RNAz (-z)
#
# # input file format:
# [1-12] sort + TPM + RNAz (at least 14 columns)
#
#########################################

use strict;
use warnings;
use Getopt::Std;

filt_sRNAs();
exit(1);

sub filt_sRNAs {
    my %opts = (l => 40, e => 20);
    getopts("l:e:", \%opts);
    usage() if (@ARGV != 1);
    my $infile = shift(@ARGV);
    my $t = read_file($infile);
    my %info = %$t;
    for my $i (sort keys %info) {
        my ($length, $tpm, $zscore) = (split /\s+/, $info{$i})[2,12,13];
# add RNAz step
#        print $info{$i} . "\n" if($zscore > 0.5);
        next if($length < $opts{l}); # length cut-off: 40 
        next if($tpm < $opts{e});    # length cut-off: 20
        print $info{$i} . "\n";        
    }
}

sub read_file {
    my $in = shift(@_);
    my %r = ();
    open my $fh_in, "< $in" or die "cannot open $in, $!\n";
    while(<$fh_in>) {
        chomp;
        next if(/^\#|^\s*$/); # blank lines
        my @tabs = split /\s+/;
        warn("[line-$.] of file $in contains less than 14 columns\n") if(@tabs < 14);
        $r{$tabs[0]} = $_;
    }
    close $fh_in;
    return (\%r);
}

sub usage{
    die(qq/
Usage: chk_filt_sRNAs.pl [options] <infile>

Options: -l <INT>   : Length cutoff. [40]
         -e <INT>   : TPM cutoff, [20]
         <infile>   : input file with at least 14 columns.

Example:
chk_filt_sRNAs.pl -l 40 -e 20 a.txt > a_filt.txt
            \n/);
}
