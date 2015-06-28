#!/usr/bin/env perl

##################################
# Classify tags from sort files
##################################

use warnings;
use strict;

sub help{
    print STDERR <<EOF;

Usage: perl sort2candi_v1.pl  infile

Note: the rules to select sRNA
1. at IGR or AS.
2. >=100 bp from 5' ORF and >=60 bp from 3' ORF.

EOF
}

my $infile = shift;
if(!defined $infile){
    &help();
    exit(1);
}

open F,$infile or die;
my ($total_input, $total_sRNA, $total_AS, $candi_AS, $total_IGR, $candi_IGR, $total_UTR, $total_IM, $total_PM, $total_other) = (0,0,0,0,0,0,0,0,0,0);
my ($print_sRNA, $print_UTR, $print_IM, $print_other) = ("", "", "", "");

while(<F>){
    $total_input ++;
    my $line = $_;
    chomp;
    my @lines = split /\t/;
    my ($gap_1, $gap_2, $dir, $des) = ($lines[7], $lines[9], $lines[10], $lines[11]);

    if($des =~ /AS|AM/i){
        $total_AS ++;
        if(($dir eq '/+/+/+/' && $gap_1 >= 100 && $gap_2 >= 60) || 
           ($dir eq '/-/-/-/' && $gap_1 >= 60 && $gap_2 >= 100) ||
           ($dir eq '/+/+/-/' && $gap_1 >= 100) ||
           ($dir eq '/-/-/+/' && $gap_1 >= 60)  ||
           ($dir eq '/-/+/+/' && $gap_2 >= 60)  ||
           ($dir eq '/+/-/-/' && $gap_2 >= 100) ||
           ($dir eq '/+/-/+/') || 
           ($dir eq '/-/+/-/')){
            $print_sRNA .= $line;
            $candi_AS ++;
        }else{
            $print_UTR .= $line;
            $total_UTR ++;
        }
    }elsif($des =~ /IGR/i){
        $total_IGR ++;
        if(($dir eq '/+/+/+/' && $gap_1 >= 100 && $gap_2 >= 60) ||
           ($dir eq '/-/-/-/' && $gap_1 >= 60 && $gap_2 >= 100) ||
           ($dir eq '/+/+/-/' && $gap_1 >= 100) ||
           ($dir eq '/-/-/+/' && $gap_1 >= 60)  ||
           ($dir eq '/-/+/+/' && $gap_2 >= 60)  ||
           ($dir eq '/+/-/-/' && $gap_2 >= 100) ||
           ($dir eq '/+/-/+/') ||
           ($dir eq '/-/+/-/')){
            $print_sRNA .= $line;
            $candi_IGR ++;
        }else{
            $print_UTR .= $line;
            $total_UTR ++;
        }    
    }elsif($des =~ /IM/i){
        $print_IM .= $line;
        $total_IM ++;
    }elsif($des =~ /PM|CM/i){
        $print_UTR .= $line;
        $total_PM ++;
    }else{
        $print_other .= $line;
        $total_other ++;
    }
}
close F;

$infile =~ /(.*)\.txt/;
my $head_name = (defined $1)?$1:"out";

open OUT,">$head_name\_sRNA\.txt" or die;   print OUT $print_sRNA;  close OUT;
open OUT,">$head_name\_UTR\.txt"  or die;   print OUT $print_UTR;   close OUT;
open OUT,">$head_name\_IM\.txt"   or die;   print OUT $print_IM;    close OUT;
open OUT,">$head_name\_OTR\.txt"  or die;   print OUT $print_other; close OUT;

$total_sRNA = $candi_AS + $candi_IGR;
my $stat = "Stat the total sRNAs:
ID\tTotal\tOutput
AS\:\t$total_AS\t$candi_AS
IGR\:\t$total_IGR\t$candi_IGR
PM\:\t$total_PM\t0
IM\:\t$total_IM\t0
Other\:\t$total_other\t0
Sum\:\t$total_input\t$total_sRNA

UTR\:\t$total_UTR\t0";

open OUT,">$head_name\.stat" or die;           
print OUT $stat . "\n";        
close OUT;
