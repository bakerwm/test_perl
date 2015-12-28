#!/bin/bash/perl -w
use warnings;
use strict;

#Convert the tab-separated file to fasta file.
#
#Input: Length Count Seq
#Output: 
#>t000001 Count
#Seq

sub usage {
    die("Usage: txt_to_fa.pl <input.txt>\n");
}

usage() if(@ARGV == 0 && -t STDIN);

my $count = 1;
while(<>) {
    chomp;
    my @tabs = split /\s+/, $_;
    if(@tabs < 3) {
        die("[format not known]\n");
    }
    if(! $tabs[1] =~ /^\d+$/ || ! $tabs[2] =~ /^[ATCGN]+$/) {
        die("[$_] not match: Len NUM Seq, format\n");
    }
    my $id = sprintf ">t%08d $tabs[1]", $count;
    print $id . "\n";
    print $tabs[2] . "\n";
    $count ++;
}



### END OF FILE ###
