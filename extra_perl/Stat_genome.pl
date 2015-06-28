#!/usr/bin/perl -w
use strict;
use warnings;
use Data::Dumper;

sub help{
    print STDERR <<EOF
Usage: perl $0  gff  match_genome.txt

Option:
    gff                 : in GFF3 format.
    match_genome.txt    : match_genome file.
EOF
}

my $gff          = shift or die &help();
my $match_genome = shift or die &help();
my $cutoff       = 0.1;

my (%mR, %rR, %tR, %Coding);
open IN,$gff or die;
while(<IN>){
    chomp;
    next unless(/\tgene\t/);
    /locus_tag=(\w+)/;
    my $id = $1;
    my @ps = split/\t/;

    ($ps[3], $ps[4]) = ("-$ps[4]", "-$ps[3]") if ($ps[6] eq "-");
# all Coding regions
    for(my $i=$ps[3]; $i<=$ps[4]; $i++){
        $Coding{$i} = 1;
    }
# mRNA
    if($id =~ /Rv\d/){
        for(my $i=$ps[3]; $i<=$ps[4]; $i++){
            $mR{$i} = 1;
        }
    }
# tRNA    
    if($id =~ /Rvnt\d/){
        for(my $k=$ps[3]; $k<=$ps[4]; $k++){
            $tR{$k} = 1;
        }
    }
# rRNA    
    if($id =~ /Rvnr\d/){
        for(my $m=$ps[3]; $m<=$ps[4]; $m++){
            $rR{$m} = 1;
        }
    }
}
close IN;

my $match_genome_uniq  = 0; 
my $match_genome_total = 0; 
my $match_mRNA_uniq    = 0; 
my $match_mRNA_total   = 0; 
my $match_tRNA_uniq    = 0; 
my $match_tRNA_total   = 0; 
my $match_rRNA_uniq    = 0; 
my $match_rRNA_total   = 0; 
my $match_AS_uniq      = 0; 
my $match_AS_total     = 0; 
my $match_IGR_uniq     = 0; 
my $match_IGR_total    = 0;
my $match_coding_uniq  = 0;
my $match_coding_total = 0;

open IN,$match_genome or die;
while(<IN>){
    my @ps = split/\t/;
    my $read_length = $ps[3] - $ps[2] + 1;
    my ($BE, $EN, $STR, $EXP) = ($ps[2], $ps[3], $ps[4], $ps[6]);
    ($BE, $EN) = ("-$EN", "-$BE") if($ps[6] eq "-");

    $match_genome_uniq  ++;
    $match_genome_total += $EXP;

    my $check_mR = 0;
    my $check_rR = 0;
    my $check_tR = 0;
    my $check_AS = 0;
    my $check_IGR= 0;
    my $over_lap = 0;

    for(my $n=$BE; $n<=$EN; $n++){
        if(exists $Coding{$n}){
            $over_lap ++;        
            $check_mR ++ if(exists $mR{$n});
            $check_rR ++ if(exists $rR{$n});
            $check_tR ++ if(exists $tR{$n});
        }elsif(exists $Coding{"-$n"}){
            $check_AS  ++;
        }else{
            $check_IGR ++;
        }
    }
    
    if($over_lap/$read_length >= $cutoff){
        if($check_mR/$read_length >= $cutoff){
            $match_mRNA_uniq  ++;
            $match_mRNA_total += $EXP;
        }elsif($check_rR/$read_length >= $cutoff){
            $match_rRNA_uniq  ++;
            $match_rRNA_total += $EXP;
        }elsif($check_tR/$read_length >= $cutoff){
            $match_tRNA_uniq  ++;
            $match_tRNA_total += $EXP;
        }else{
            $match_coding_uniq  ++;
            $match_coding_total += $EXP;
        }
    }elsif($check_AS/$read_length >= $cutoff){
        $match_AS_uniq  ++;
        $match_AS_total += $EXP;
#    }elsif($check_IGR/$read_length >= $cutoff){
    }else{
        $match_IGR_uniq  ++;
        $match_IGR_total += $EXP;
    }
}

print "Genome\t$match_genome_uniq\t$match_genome_total
mRNA\t$match_mRNA_uniq\t$match_mRNA_total
tRNA\t$match_tRNA_uniq\t$match_tRNA_total
rRNA\t$match_rRNA_uniq\t$match_rRNA_total
AS\t$match_AS_uniq\t$match_AS_total
IGR\t$match_IGR_uniq\t$match_IGR_total
Other\t$match_coding_uniq\t$match_coding_total
";
