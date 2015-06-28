#!/usr/bin/perl  -w
use strict;
use warnings;
use Getopt::Std;

#################################################
# Stat length distribution for FASTA/Sorted file.
# Date: 2012-12-01
#################################################

sub help{
    print STDERR <<EOF
Usage: perl  length_dis.pl  -i  clean.fa  -o  length_dis.txt  -s  length.stat

Note:
    i       : Input file. FASTA/FASTQ/Sorted file.
    o       : Output result.
    s       : State of length distribution.
EOF
}

my %options = ();
getopts("i:o:s:", \%options);
if(!defined $options{i}){
    &help();
    exit(1);
}

$options{o} = (defined $options{o})?$options{o}:'length_dis.txt';
$options{s} = (defined $options{s})?$options{s}:'length.stat';

# Check the file format
open F,$options{i} or die;
my $format;
while(<F>){
    if(/^>/){
        $format = "FASTA";
        last;
    }elsif(/^@/){
        $format = "FASTQ";
        last;
    }else{
        my @t = split(/\t/);
        if(@t >= 6 && $t[3] =~ /^\d+$/ && $t[4] =~ /^\d+$/ && $t[5] =~ /\+|\-/){
            $format = "Sorted";
            last;
        }else{
            &error();
            exit;
        }
    }
}
close F;

open F, $options{i} or die;
open OUT,"> $options{o}" or die;
open STA,"> $options{s}" or die;
my %fa;
# Cal fasta file.
if($format eq 'FASTA'){
    $/ = ">";
    while(<F>){
        $_ =~ s/>//g;
        my @lines = split(/\n/, $_);
        next if(@lines <= 1);
        shift(@lines);  # read id
        my $seq = join"",@lines;
        my $length = length($seq);
        $fa{$length} ++;
    }
    $/ = "\n";
# Cal FastQ file.
}elsif($format eq 'FASTQ'){
    $/ = "@";
    while(<F>){
        $_ =~ s/@//g;
        my @lines = split(/\n/, $_);
        next if(@lines <= 1);
        my $seq = $lines[1];
        my $length = length($seq);
        $fa{$length} ++;
    }
    $/ = "\n";
# Cal Sorted file.
}elsif($format eq 'Sorted'){
    foreach my $m (<F>){
        my @t = split(/\t/, $m);
        next if($t[3] =~ /\D/ || $t[4] =~ /\D/);
        my $length = $t[4] - $t[3] + 1;
        $fa{$length} ++;
    }
}else{
    &error();
    exit;
}
close F;

# print out file.
my ($total, $mean, $num) = (0,0,0);
foreach my $n (sort{$a<=>$b} keys %fa){
    $num += $fa{$n};
    $total += $n * $fa{$n};
    print OUT $n,"\t",$fa{$n},"\n";
}

my ($min, $max) = &min_max(keys %fa);
$mean = sprintf "%.2f",($total/$num);
my $stat = "Total Reads: $num
Maximum length: $max
Minimum length: $min
Mean length: $mean";
print STA $stat,"\n";
close OUT;
close STA;

sub error{
    print STDERR <<EOF
Input file should be in FASTA or Sorted format.

Sorted format: [at least 6 lines]
<id*> <Exp> <Length> <Begin*> <End*> <Strand*> <...>
EOF
}

sub min_max{
    my $min = my $max = shift;
    foreach(@_){
        $min = ($min < $_)?$min:$_;
        $max = ($max > $_)?$max:$_;
    }
    return ($min, $max);
}
