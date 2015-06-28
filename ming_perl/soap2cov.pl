#!/usr/bin/perl -w
use warnings;
use strict;

sub help{      
    print STDERR  <<EOF
Usage: perl  $0  genome.fa  Rv.coverage  Rv-1.soap  Rv-2.soap ...

Note:
    genome.fa       : The reference genome in fasta format.
    Rv              : The output filename.
    Rv-1.soap       : The soap SE or soap PE after merged,
                      multiple soap files are acceptable.

This script is designed to calculate the coverage for each base
using the soap files.

EOF
}

my $fa   = shift or die &help(); # genome fasta
my $name = shift or die &help(); # name of output file.
my @soap = @ARGV or die &help(); # soap files

open F,$fa or die;
my ($genome_len, $genome_id);
while(<F>){
    chomp;
    if(/\>/){
        my @id = split/\s+/;
        $id[0] =~ s/\>//;
        $genome_id = $id[0];
        next;
    }
    my $length = length($_);
    $genome_len += $length;
}
close F;

my %hash;
foreach my $file (@soap){
    open F, $file or die;   # SOAP result, total TO uniq;
    while(<F>){
        chomp;
        my @ps = split/\t/;
        my $read = join"\,",($ps[5], $ps[6], $ps[8]);   # len   +   position
        $hash{$read} ++;
    }
    close F;
}

my %cov;
foreach my $i(keys %hash){  # Cal coverage
    my @ps = split/\,/,$i;
    for(my $k=$ps[2]; $k<=($ps[2]+$ps[0]-1); $k++){
        $cov{$ps[1]}->{$k} += $hash{$i};
    }
}

my ($out_n, $out_p);
for(my $j=1; $j<=$genome_len; $j++){
    my $cov_n = (exists $cov{'-'}->{$j})?$cov{'-'}->{$j}:'0';
    my $cov_p = (exists $cov{'+'}->{$j})?$cov{'+'}->{$j}:'0';

    $out_n .= "$genome_id\t$j\t$cov_n\n";
    $out_p .= "$genome_id\t$j\t$cov_p\n";
}

open O1, "> $name\.n" or die;
open O2, "> $name\.p" or die;
print O1 $out_n;
print O2 $out_p;
close O1;
close O2;
