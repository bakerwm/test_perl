#!/usr/bin/env perl

##################################################
# Rename the ncRNAs according to the publication:
#
# Lamichhane G, et al., Definition and annotation 
# of (myco)bacterial non-coding RNA, Tuberculosis 
# (2012), 
# http://dx.doi.org/10.1016/j.tube.2012.11.010
##################################################

use strict;
use warnings;

use Data::Dumper;

my $id_prefix = 'ncRv';
my $length_prefix = 2;
#my $id_prefix = 'ncMRA';
#my $length_prefix = 4;

my $usage="perl  $0  IN.ptt  list.txt  >out.txt\n";
my $ptt=shift or die "Input ptt file:\n$usage";
my $infile=shift or die "Input txt file:\n$usage";

open F, $ptt or die "Cannot open file $ptt: $!";
my %Gene;
my $lab=1;
while(<F>){
    chomp;
    next unless(/^\d+\.\.\d+/);
    my @ps=split/\t/;
    $Gene{$lab}=$ps[5];
    $lab++;
}
close F;

my %ID;
foreach my $i(sort{$a<=>$b} keys %Gene){
    my $pre=$i-1;
    my $nxt=$i+1;
    my $preGene=(exists $Gene{$pre})?$Gene{$pre}:$Gene{$i};
    my $nxtGene=(exists $Gene{$nxt})?$Gene{$nxt}:$Gene{$i};
    $ID{$Gene{$i}}="$preGene\t$nxtGene";
}

open F, $infile or die "Cannot open file $infile: $!";
my %IGR;
my %AS;
while(<F>){
    chomp;
    my $line=$_;
    my @ps=split/\t/,$_;
    if($ps[11] eq "IGR"){
        my $tag="$ps[6]\:$ps[8]";
        push @{$IGR{$tag}},$_;    
    }elsif($ps[11]=~/AS|PM|IM/){
        my $tag="";
        if($ps[6] eq $ps[8]){
            $tag=$ps[6];
        }else{
            if($ps[7]<=0 && $ps[9]<=0){
                $tag=$ps[6];
            }elsif($ps[7]<=0 && $ps[9]>=0){
                $tag=$ps[6];
            }elsif($ps[7]>=0 && $ps[9]<=0){
                $tag=$ps[8];
            }elsif($ps[7]>=0 && $ps[9]>=0){
                $tag=$ps[8];
            }
        }
        push @{$AS{$tag}},$_;
    }else{
        print $_,"\n";
    }
}
close F;

my @fix=qw/A B C D E F G H I J K L M N O P Q R S T U V W X Y Z AA AB AC AD AE AF AG AH AI AJ AK AL AM AN AO AP AQ AR AS AT AU AV AW AX AY AZ/;
foreach my $g(sort keys %IGR){
    my $rank=0;
    foreach my $lib(@{$IGR{$g}}){
        my @ma=split/\t/,$lib,13;
        my $str=$ma[5];
        my ($R1,$R2)=split/\:/,$g;
#        my $head="ncRa1";       
#        my $head="ncRv1";
        my $head = $id_prefix . '1';
        my $num=&RvNum($R1);
        my $ord=$fix[$rank];
        my $tail=($str eq "-")?"c":"";
        my $NewID=$head.$num.$ord.$tail;
        ($ma[12],$ma[0])=($ma[0],$NewID);
        my $p_out=join"\t",@ma;
        print $p_out,"\n";

        $rank++;
    }
}

foreach my $a(sort keys %AS){
    my $rank=0;
    foreach my $lin(@{$AS{$a}}){
        my @mb=split/\t/,$lin,13;
        my $str=$mb[5];
#        my $head="ncRa";
#        my $head="ncRv";
        my $head = $id_prefix;
        my $num=&RvNum($a);
        my $ord=$fix[$rank];
        my $tail=($str eq "-")?"c":"";
        my $NewID=$head.$num.$ord.$tail;
        ($mb[0],$mb[12])=($NewID,$mb[0]);
        my $p_out=join"\t",@mb;
        print $p_out,"\n";
        
        $rank++; 
    }
}


# For H37Ra
sub RvNum{
    my $in=shift;
#    $in=~s/[a-zA-Z]$//;
    $in =~ s/c$//;
    my $end=length($in) - $length_prefix;
    my $num=substr($in, $length_prefix ,$end);
    return $num;
}

sub usage {
    die("
Usage: Name_ncRNA.pl [options] <in.txt>

Options: -a <STR>   : annotation file (PTT)
         -p <STR>   : the prefix of the new name

         <in.txt>   : tab-separated file

Notes:
1. new name would be in the following formats:
ncRv1234A ncRv1234Ac ncRv11234
\n");
}

