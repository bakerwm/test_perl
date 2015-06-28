#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;

########################################
# Readme:
# The file list.pos.out consist of the following columns:
# col1-3:  Geneid, cov, length;
# col4-6:  The left edge of sRNA, the right edge sRNA. 
#          To determin frame of the figure.
# col7-10: The minimum position of matched sRNA, the maximum 
#          position of matched sRNAs,(exclude mRNA/sRNA)
# col11- : The coordination for each mapped sRNA.
# col last two: The 5' and 3' gap.
##########################################

sub help{
    print STDERR <<EOF

Usage: perl  match.pl  -a total_sRNA.txt  -i input.txt  -c 0.5 -o match_input.txt
Options:
    -a <file>       : The total sRNA file.
    -i <file>       : The input sequence file.
    -o <file>       : Ouput file.
    -c [0-1]        : The cutoff for overlap with each sequence.
                      default is 0.5.
    -h              : Show this help.

Note:
1. Both input file should be in this format:
ID* Exp Length  Begin*  End*    Strand* (* are required)
2. Cutoff is no more than 1.

Description:
1. This script is designed for detecting sRNAs that have overlap
   with input sequence according to genome coordination.
2. The overlap between two sequences is no less than  the cutoff 
   for each sequence.
        
EOF
}

my ($total_sRNA, $infile, $cutoff, $outfile, $help);
&GetOptions(
            "a=s"   =>  \$total_sRNA,
            "i=s"   =>  \$infile,
            "c=s"   =>  \$cutoff,
            "o=s"   =>  \$outfile,
            "h"     =>  \$help,
           );

$cutoff = 0.5 unless(defined $cutoff);

if($help or not $total_sRNA or not $infile or not $outfile){
    &help();
    exit(1);
}

my (%total);
open (F,"<$total_sRNA") or die;   # For total sRNA.
while(<F>){          
    chomp;
    my @ps = split (/\t/);
    @{$total{$ps[0]}} = @ps;
}
close F;

# Search the sRNAs have overlap with input sequence.
my (%input, %hit);
open (F, "<$infile") or die; # Input sequence
while(<F>){
    chomp;
    my @ps = split (/\t/);
    @{$input{$ps[0]}} = @ps; #
    my ($In_begin, $In_end, $In_strand) = ($ps[3], $ps[4], $ps[5]);
    my $In_length = $In_end - $In_begin + 1;
    foreach my $id(sort keys %total){
        my @t=@{$total{$id}};
        my ($sRNA_begin, $sRNA_end, $sRNA_strand) = ($t[3], $t[4], $t[5]);
        my $sRNA_length = $sRNA_end - $sRNA_begin + 1;
        next unless($sRNA_strand eq $In_strand);

        my $match  = my $gap = 0;
        if($sRNA_begin < $In_begin && $sRNA_end > $In_begin){
            $gap   = $sRNA_end - $In_begin + 1;
            $match = ($gap > $In_length)?$In_length:$gap;
        }elsif($sRNA_begin >= $In_begin && $sRNA_begin < $In_end){
            $gap   = $In_end - $sRNA_begin + 1;
            $match = ($gap > $sRNA_length)?$sRNA_length:$gap;
        }else{
            next;
        }
# Make sure overlap larger than 1/2 of each fragement.
        next unless($match/$In_length >= $cutoff || $match/$sRNA_length >= $cutoff || $match >= 20);
        push(@{$hit{$ps[0]}}, $t[0]);
    }
}
close F;

#############################################
# Count the max number of sRNAs in one line.
#############################################
my $max_sRNA = 0;
foreach my $g (sort keys %hit){
    $max_sRNA = ($max_sRNA > @{$hit{$g}})?$max_sRNA:@{$hit{$g}}; 
}

#############################################
# Determin the edges for sRNAs match one 
# sequence.
#############################################
my $outfile2 = $outfile;
$outfile2 =~ s/\.txt/\.list/g;

open (O1,"> $outfile" ) or die;
open (O2,"> $outfile2") or die;

foreach my $g (sort keys %input){  #  
    next unless(exists $hit{$g});  # skip the sequences that do not match any sRNAs.
    my @b = @{$input{$g}};    # The coordination of mRNA/sRNA: "Rv0001, exp, 1524, 1, 1524, +"
    my ($In_begin, $In_end, $In_strand) = ($b[3], $b[4], $b[5]);
#    $m_b=$b[3]; $m_e=$b[4]; $m_s=$b[5];  
    
    my @t = @{$hit{$g}};  # contain the sRNAs overlap with mRNA/sRNA.
    my $min_edge = $total{$t[0]}[3];    # t[0] is the first sRNA match input sequence.
    my $max_edge = $total{$t[0]}[4]; 
    
    my @pos = ();
    for(my $i = 0; $i < $max_sRNA; $i++){  # keep each row in the same length.
# add position of input sequence to keep all rows in the same length.
        my @sRNA_match = (defined ($t[$i]))?@{$total{$t[$i]}}:@b;  
        push @pos, ($sRNA_match[3], $sRNA_match[4], $sRNA_match[5]);
        next if($sRNA_match[0] eq $g);    # if this is the added sRNA, next.
        ($min_edge, $max_edge) = &min_max($min_edge, $max_edge, $sRNA_match[3], $sRNA_match[4]);
    }

    my $gap_5end = $In_begin - $min_edge; 
    my $gap_3end = $max_edge - $In_end;
    ($gap_5end, $gap_3end) = ($gap_3end, $gap_5end) if($In_strand eq "-");

    my ($figure_left, $figure_right) = &min_max($min_edge, $max_edge, $In_begin, $In_end);
    
    my $line_1 = join"\t",($g, "exp", "len", $figure_left, $figure_right, $In_strand, $min_edge, $max_edge, $In_strand, $In_begin, $In_end, $In_strand, @pos, $gap_5end, $gap_3end);
    my $line_2 = join"\t",($g, @{$hit{$g}});

    print O1 $line_1,"\n";
    print O2 $line_2,"\n";
        
    @pos = ();
}

my $max_column = ($max_sRNA + 4)*3 + 2;
print "The most sRNAs and lines are:\t$max_sRNA\t$max_column\n";
#print "The most sRNA matched in one line:\t$max_sRNA\nThe number of columns in pos.txt file:\t$max_column\n";

close O1;
close O2;

sub min_max{
    my $a = my $b = shift(@_);
    foreach(@_){
        $a=($a>$_)?$_:$a;
        $b=($b<$_)?$_:$b;
    }
    return($a, $b);
}
