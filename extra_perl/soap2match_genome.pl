#!/usr/bin/perl  -w
use warnings;
use strict;
use Getopt::Std;

sub help{
    print STDERR <<EOF
Usage: perl  soap2match_genome.pl  -i  in.soap  -o match_genome.txt

Noete:
    i       : SE soap file or PE merged file.
    o       : match_genome file

EOF
}

my %options = ();
getopts("i:o:", \%options);
if(!defined $options{i} || !defined $options{o}){
    &help();
    exit(1);
}

my %seq;
open F, $options{i} or die;
foreach my $i(<F>){
    my @tabs = split(/\t/,$i);
    $seq{$tabs[1]}->{num} ++;
    $seq{$tabs[1]}->{info} = join"\t",($tabs[5], $tabs[6], $tabs[7], $tabs[8]);  # length, str, ref, begin.
}
close F;

open OUT,"> $options{o}" or die;
my $num = 1;
my $print_out;
foreach my $k (keys %seq){
    my @info = split(/\t/, $seq{$k}->{info});
    my $end  = $info[0] + $info[3] - 1;
    my $tag  = sprintf "%08d",$num;
    my $out  = join"\t",("t$tag", $info[2], $info[3], $end, $info[1], $k, $seq{$k}->{num});
#    print OUT $out,"\n";
    $print_out .= "$out\n";
    $num ++;
}
print OUT $print_out;

close OUT;
