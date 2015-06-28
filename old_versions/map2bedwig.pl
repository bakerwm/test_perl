#! /usr/bin/perl -w
# Filename: map2bedwig.pl
# Convert bowtie map file into *.bed and *.wig file
# Note: reads are in fasta format as follow:
# >seqid1_x56
# CAGATCGATCGATCGAGCGCGACT
# >seqid2_x48
# GTCGATGAGCGACCGCGCGCGCGG
# Here 56 and 48 after x are the reads numbers.
# qinhu.wang (at) gmail.com

use strict;

my $usage = "  
Usage: $0
MAP: bowtie map file without ¨Csuppress
SEP: ID Seperator of reads_id and reads_num\n"; 

print $usage;

my $map = $ARGV[0] || "bowtie.map";
my $sep = $ARGV[1] || "_x";
my %bed = ();
my %wig = ();

open (MAP, $map) or die "Cannot open file $map: $!\n";
while (<MAP>) {
my ($index, $strand, $refseqid, $start, $read, $discard) = split /\t/, $_, 6;
my ($idx, $num) = split /$sep/, $index;
my $end = $start + length($read);
$bed{$index} = $refseqid."\t".$start."\t".$end."\t".$strand;
$wig{$index} = $refseqid."\t".$start."\t".$end."\t".$num;
}
close MAP; 

open (BED, ">map.bed") or die "Cannot create file map.bed: $!\n";
open (WIG, ">map.wig") or die "Cannot create file map.wig: $!\n";
foreach (sort keys %bed) {
print BED $bed{$_}, "\n";
print WIG $wig{$_}, "\n";
}
close BED;
close WIG;