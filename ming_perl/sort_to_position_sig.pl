#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Std;
use Data::Dumper;

######################################################################
# This script is designed to determin the genome position of input   #
# sRNAs.                                                             #
# Date: 2012-12-22                                                   #
#                                                                    #
# update: 2013-12-21                                                 #
#   Line-82                                                          #
#   Change the regular expression for coverage number in case they   #
#   are in scientific expression (1e+2).                             #
#   replace [\w+] by [.*].                                           #
######################################################################

my %opts = (t => "gene");
getopts("f:t:g:", \%opts);
die "Usage: perl $0 [-g] <genome.gff> [-f] <ref.fa> input.txt > out.txt" if (@ARGV != 1);

my $in_gff = $opts{g};
my $in_fa  = $opts{f};
my $infile = shift;
################################################################################
# Readin exclude list
my %ex_list = (); # &read_txt($exclude) if (defined $exclude);

# Read fa length;
my $genome_length = 0;
open F, "< $in_fa" or die "$!";
while(<F>){
    chomp;
    next if(/\^>/);
    $genome_length += length($_);
}
close F;

# Read gff file. 
open F, $in_gff or die "Cannot open $in_gff $!\n";
my @lines = <F>;
close F;

my %gff = ();
## Parse the gff file
foreach (@lines){
    next if(/^\#/);
    my($g_begin,$g_end,$g_strand)=(split/\t/)[3,4,6];
    my $g_name;
    if(/\t$opts{t}\t/) {
#    if(/\tgene\t/){
        if(($g_name) = $_ =~ /locus_tag=(\w+)/){
        }elsif(($g_name) = $_ =~ /Name=(\w+)/){
        }elsif(($g_name) = $_ =~ /ID=gene\:(\w+)/){
        }elsif(($g_name) = $_ =~ /GeneID\:(\d+)/){
        }else{
        }
#        if(/^K_051809/){ ## For B42.gff Only
#            ($g_name) = /Name=(\w+)/;
#        }elsif(/Name=\w+/){
#            ($g_name) = /locus_tag=(\w+)/ if(/locus_tag=/);  ### For locus tag
#        }else{
#            ($g_name) = /locus_tag=(\w+)/; ## NCBI gff version 1.14
#        }
        $gff{$g_begin} = join "\t", ($g_name,$g_begin,$g_end,$g_strand);
    }elsif(/\tncRNA\w+\t/){
        if(($g_name) = $_ =~ /locus_tag=(\w+)/){
        }elsif(($g_name) = $_ =~ /Name=(\w+)/){
        }elsif(($g_name) = $_ =~ /ID=gene\:(\w+)/){
        }elsif(($g_name) = $_ =~ /GeneID\:(\d+)/){
        }else{
        }
        $ex_list{$g_begin} = join "\t", ($g_name,$g_begin,$g_end,$g_strand); 
    }
}

# Delete known ncRNAs
foreach my $e (keys %ex_list) {
#    delete $gff{$e};
}

# Sort genes according to BEGIN position.
my $num = 1;
my %GENE = ();
foreach my $i(sort{$a<=>$b} keys %gff){
    $gff{$i} .= "\t$num";
    $GENE{$num} = $gff{$i};
    $num ++;
}
my $total_gene_number = (keys %gff);

# Read sRNA candidate file.
open F, $infile or die;
my @InLists = <F>;
close F;

#open OUT,"> $opts{o}" or die;
foreach my $j (@InLists){
    chomp ($j);
    $j =~ s/\n|\r//;
    my ($ca_id, $ca_chr, $ca_begin_pos, $ca_end_pos, $ca_str) = (split /\t/, $j)[0,1,3,4,5];
    my $ca_length = $ca_end_pos - $ca_begin_pos + 1;
#   my ($ca_id,$ca_str,$ca_begin_pos,$ca_max_cov,$ca_end_pos,$ca_length) = $j=~/(^\w+)\t([+,-])\t(\d+)\:.*\t\w+\:(.*)\t(\d+)\:.*\t(\d+)/; # Rv0001 + 1:Cov Max:1000 1524:Cov 1524
    my $ca_info = join"\t",($ca_id, $ca_chr, $ca_length, $ca_begin_pos, $ca_end_pos, $ca_str);
# Search the neighbor genes.
    my ($pre_gene, $next_gene) = ('Start','End');
    my ($pre_gene_order, $next_gene_order) = (1, $total_gene_number);
    my ($pre_gene_begin, $pre_gene_end, $pre_gene_str)    = (1, 1, '+');
    my ($next_gene_begin, $next_gene_end, $next_gene_str) = ($genome_length, $genome_length, '+');
    # Find previous gene.
    for(my $k=$ca_begin_pos; $k>=1; $k--){
        if(exists $gff{$k}){
            ($pre_gene, $pre_gene_begin, $pre_gene_end, $pre_gene_str, $pre_gene_order) = split(/\t/, $gff{$k});
            last;
        }else{
            next;
        }
    }
    # Find next gene.
    for(my $k=$ca_begin_pos; $k<=$genome_length; $k++){
        if(exists $gff{$k}){
            ($next_gene, $next_gene_begin, $next_gene_end, $next_gene_str, $next_gene_order) = split(/\t/, $gff{$k});
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
    print $ca_info . "\t" . $out . "\n";
}

### Subroutines ###
sub read_txt{
    my $in = shift(@_);
    my %list = ();
    open my $fh, "< $in" or die "$!";
    while(<$fh>) {
        chomp;
        my $id = (split /\s+/)[0];
        $list{$id} = $_;
    }
    return %list;
}

