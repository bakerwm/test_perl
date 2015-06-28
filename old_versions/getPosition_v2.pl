#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Std;
use Data::Dumper;

######################################################################
# This script is designed to determin the genome position of input   #
# sRNAs.                                                             #
# Date: 2012-12-22
# update: 2013-12-21
######################################################################

sub help{
    print STDERR <<EOF
Usage: perl  getPosition.pl  -i  H37Rv_sRNA.temp  -f  NC_000962.gff  -o H37Rv_sRNA.txt

Note:
    i       : Input file is tab-delimited with 6 lines.
              <id> <strand> <begin:cov> <max:cov> <end:cov> <length>
    f       : The gff annotation file.
    o       : The output file.
EOF
}

my %opts = ();
getopt("i:f:o:", \%opts);
if(!defined $opts{i} || !defined $opts{f} || !defined $opts{o}){
    &help();
    exit(1);
}

# Read GFF file.
open F, $opts{f} or die;
my %GFF;
my $genome_length;
foreach my $i (<F>){
    if($i =~ /^##sequence-region/){
        my @lines = split(/\s+/, $i);
        $genome_length = $lines[3];
    }
    next unless($i =~ /\tgene\t/);
    my @t = split(/\t/, $i);
    my $gene_tag;
    if($opts{f} =~ /B42|B36|BT/i){
        $i =~ /Name=(\w+)/;
        $gene_tag = $1;
    }else{
        $i =~ /locus_tag=(\w+)/;
        $gene_tag = $1;
    }
    $GFF{$t[3]} = join"\t",($gene_tag, $t[3], $t[4], $t[6]);
}
close F;

# Sort genes according to BEGIN position.
my $num = 1;
my %GENE = ();
foreach my $i(sort{$a<=>$b} keys %GFF){
    $GFF{$i} .= "\t$num";
    $GENE{$num} = $GFF{$i};
    $num ++;
}
my $total_gene_number = (keys %GFF);

# Read sRNA candidate file.
open F, $opts{i} or die;
open OUT,"> $opts{o}" or die;
foreach my $j (<F>){
    $j =~ /(^\w+)\t(\+|\-)\t(\d+)\:\w+\t\w+\:(\w+)\t(\d+)\:\w+\t(\d+)/; # Rv0001 + 1:Cov max:Cov 1524:Cov 1524
    my ($ca_id, $ca_str, $ca_begin_pos, $ca_max_cov, $ca_end_pos, $ca_length) = ($1, $2, $3, $4, $5, $6);
    my $ca_info = join"\t",($ca_id, $ca_max_cov, $ca_length, $ca_begin_pos, $ca_end_pos, $ca_str);
# Search the neighbor genes.
    my ($pre_gene, $next_gene) = ('Start','End');
    my ($pre_gene_order, $next_gene_order) = (1, $total_gene_number);
    my ($pre_gene_begin, $pre_gene_end, $pre_gene_str)    = (1, 1, '+');
    my ($next_gene_begin, $next_gene_end, $next_gene_str) = ($genome_length, $genome_length, '+');
    # Find previous gene.
    for(my $k=$ca_begin_pos; $k>=1; $k--){
        if(exists $GFF{$k}){
            ($pre_gene, $pre_gene_begin, $pre_gene_end, $pre_gene_str, $pre_gene_order) = split(/\t/, $GFF{$k});
            last;
        }else{
            next;
        }
    }
    # Find next gene.
    for(my $k=$ca_begin_pos; $k<=$genome_length; $k++){
        if(exists $GFF{$k}){
            ($next_gene, $next_gene_begin, $next_gene_end, $next_gene_str, $next_gene_order) = split(/\t/, $GFF{$k});
            last;
        }else{
            next;
        }
    }
    my ($L1, $L2, $L3, $L4, $R1, $R2, $R3, $R4, $R7, $R8);
    $L1 = $ca_begin_pos - $pre_gene_begin;
    $L2 = $ca_begin_pos - $pre_gene_end;
    $L4 = $ca_end_pos - $pre_gene_end;
    $R1 = $ca_begin_pos - $next_gene_begin;
    $R3 = $ca_end_pos - $next_gene_begin;
    $R4 = $ca_end_pos - $next_gene_end;
    my($next_gene_2, $next_gene_begin_2, $next_gene_end_2, $next_gene_str_2) = ('End', $genome_length, $genome_length, '+');
    if($next_gene_order < $total_gene_number){
        my $next_gene_order_2 = $next_gene_order + 1;
        ($next_gene_2, $next_gene_begin_2, $next_gene_end_2, $next_gene_str_2) = split(/\t/, $GENE{$next_gene_order_2});
    }
    $R7 = $ca_end_pos - $next_gene_begin_2;
    $R8 = $ca_end_pos - $next_gene_end_2;
# Determin the gap/direction/description between candidate and pre/next gene.
    my ($gap_1, $gap_2, $direction, $des);
    if($L1>=0 && $L2<=0){
        if($L4<=0){
            $gap_1 = $L1; $gap_2 = -$L4;
            $next_gene = $pre_gene;
            $direction = "\/$pre_gene_str\/$ca_str\/$pre_gene_str\/";
            $des = ($ca_str eq $pre_gene_str)?'IM':'AS';
        }elsif($L4>0 && $R3<0){
            $gap_1 = $L2; $gap_2 = -$R3;
            $direction = "\/$pre_gene_str\/$ca_str\/$next_gene_str\/";
            $des = ($ca_str eq $pre_gene_str)?'PM':'AS';
        }elsif($R3>=0 && $R4<=0){
            $gap_1 = $L2; $gap_2 = -$R3;
            $direction = "\/$pre_gene_str\/$ca_str\/$next_gene_str\/";
            $des = ($ca_str ne $pre_gene_str && $ca_str ne $next_gene_str)?'AS2':'PM2';
        }elsif($R4>0 && $R7<0){
            $next_gene = $next_gene_2;
            $gap_1 = $L2; $gap_2 = -$R7;
            $direction = "\/$pre_gene_str\/$ca_str\/$next_gene_str_2\/";
            $des = ($ca_str ne $pre_gene_str && $ca_str ne $next_gene_str)?'AS2':'PM2';
        }elsif($R7>=0 && $R8<=0){
            $next_gene = $next_gene_2;
            $gap_1 = $L2; $gap_2 = -$R7;
            $direction = "\/$pre_gene_str\/$ca_str\/$next_gene_str_2\/";
            $des = ($ca_str ne $pre_gene_str && $ca_str ne $next_gene_str && $ca_str ne $next_gene_str_2)?'AS3':'PM3';
        }else{
            $next_gene = $next_gene_2;
            ($gap_1, $gap_2) = (1, 1);
            $direction = '/+/+/+/';
            $des = 'Null';
        }
    }else{
        if($R3<0){
            $gap_1 = $L2; $gap_2 = -$R3;
            $direction = "\/$pre_gene_str\/$ca_str\/$next_gene_str\/";
            $des = 'IGR';
        }elsif($R3>=0 && $R4<=0){
            $gap_1 = $L2; $gap_2 = -$R3;
            $direction = "\/$pre_gene_str\/$ca_str\/$next_gene_str\/";
            $des = ($ca_str eq $next_gene_str)?'PM':'AS';
        }elsif($R4>0 && $R7<0){
            $next_gene = $next_gene_2;
            $gap_1 = $L2; $gap_2 = -$R7;
            $direction = "\/$pre_gene_str\/$ca_str\/$next_gene_str_2\/";
            $des = ($ca_str eq $next_gene_str)?'PM2':'AS2';
        }elsif($R7>=0 && $R8<=0){
            $next_gene = $next_gene_2;
            $gap_1 = $L2; $gap_2 = -$R7;
            $direction = "\/$pre_gene_str\/$ca_str\/$next_gene_str_2\/";
            $des = ($ca_str ne $next_gene_str && $ca_str ne $next_gene_str_2)?'AS2':'PM2';
        }else{
            $next_gene = $next_gene_2;
            ($gap_1, $gap_2) = (1, 1);
            $direction = '/+/+/+/';
            $des = 'Null';
        }
    }
    my $out = join"\t",($pre_gene, $gap_1, $next_gene, $gap_2, $direction, $des);
    print OUT $ca_info,"\t",$out,"\n";
}
close F;
close OUT;
