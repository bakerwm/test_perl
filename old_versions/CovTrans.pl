#!/usr/bin/perl -w

use strict;

##############################################
# This script is designed to transform output 
# by soap.coverage -depthsingle [text, fasta]
# to tab-delimited format.
##############################################

my $usage = "Usage: perl  $0  in.coverage result.out  Genome_fna coverage
Note:
1. Infile is the result from soap.coverage -depthsingle
2. Length is the result from soap.coverage -o
3. The genome fasta file used for soap
4. The result file\n";


my $infile  = shift or die $usage;
my $length  = shift or die $usage;
my $genome  = shift or die $usage;
my $outfile = shift or die $usage;

open F,"<$length" or die;
my $genome_length = 0;
while(<F>){
    chomp;
    next unless(/^Total/i);
    /Total\:(\d+)/i;
    $genome_length = $1;
}
close F;

open F,"<$genome" or die;
my $genome_id;
while(<F>){
    chomp;
    next unless(/\>/);
    s/\>//g;
    my @ps = split(/\s/);
    $genome_id = $ps[0];
}
close F;

open F,"<$infile" or die;
my $pos   = 1;
my $p_out = "";
while(<F>){
    chomp;
    next if(/\>/);
    my @ps = split(/\s/);

    foreach(@ps){
        $p_out .= $genome_id."\t".$pos."\t".$_."\n";
        $pos ++;
    }
}
close F;

$pos --; 

if($genome_length ne $pos){
    print "Check: the genome length in $infile \($pos\) and $length \($genome_length\) does not match.\n";
    exit(1);
}

open OUT,">$outfile" or die;
print OUT $p_out;
close OUT;
