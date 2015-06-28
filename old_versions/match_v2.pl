#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Std;

################################################################
# Readme:                                                      #
# The file list.pos.out consist of the following columns:      #
# col1-3:   Geneid, cov, length;                               #
# col4-6:   The edges of input seq and hit RNAs.               #
#           Determin the frame of the figure.                  #
# col7-9:   The edges of hit RNAs.                             #
# col10-12: The coordination for input sequence.               #
# col13-..: The coordination for hit RNAs.                     #
# col last two: The 5' and 3' gap.                             #
################################################################

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
1. Both files should be in this format:
ID* Exp Length  Begin*  End*    Strand* (* are required)
2. Cutoff is no more than 1.

Description:
1. This script is designed for detecting sRNAs that have overlap
   with input sequence according to genome coordination.
2. The overlap between two sequences is no less than  the cutoff 
   for each sequence.
        
EOF
}

my %opt=();
getopt("a:i:c:o:h:",\%opt);
$opt{c}=($opt{c})?$opt{c}:"0.5";
$opt{o}=($opt{o})?$opt{o}:"out.txt";

if($opt{h} || !$opt{a} || !$opt{i}){
    &help();
    exit(1);
}

my (%total);
my $R_max=0; # The longest sequence in total file.
open (F,"$opt{a}") or die "Cannot open file $opt{a}: $!";
while(<F>){
    chomp;
    my ($beg,$end) = (split (/\t/))[3,4]; # Begin position.
    $R_max=(($end-$beg+1)>$R_max)?($end-$beg+1):$R_max;
    push @{$total{$beg}},$_;
}
close F;

my $idfile=$opt{o}.".list";
open ID,">$idfile" or die "Cannot open file $idfile: $!";
my %Cakes; # list sequences hit total files.
my $maxCount=1; # the number of sequences hit in one list seq.
open (F, "$opt{i}") or die;
while(<F>){
    chomp;
    my ($q_id,$q_beg,$q_end,$q_str)=(split/\t/)[0,3,4,5];
    my $q_len=$q_end-$q_beg+1;
    my $edge_L=$q_beg-$R_max;
    my ($FigL,$FigR)=($q_beg,$q_end);
    my (@RNAPos, @Rhits, @Rids);
    my $count=0;
    for(my $i=$q_end;$i>=$edge_L;$i--){
        next unless(exists $total{$i});
        my @R_cdt=@{$total{$i}}; # RNA candidates
        foreach my $k(@R_cdt){
            my ($s_id,$s_beg,$s_end,$s_str)=(split/\t/,$k)[0,3,4,5];
            my $s_len=$s_end-$s_beg+1;
            next unless($q_str eq $s_str); # same strand
            next unless($q_beg < $s_end && $q_end>$s_beg); # have overlap
            my $gap=&overlap($q_beg,$q_end,$s_beg,$s_end);
            my $q_pct=sprintf"%.2f",($gap/$q_len);
            my $s_pct=sprintf"%.2f",($gap/$s_len);
            next if($q_pct<$opt{c} && $s_pct<$opt{c});
            my $tmp=join"\t",($s_beg,$s_end,$s_str);
            push @Rhits,$tmp;
            push @Rids,$s_id;
            ($FigL,$FigR)=&MinMax($q_beg,$q_end,$s_beg,$s_end);
            push @RNAPos,($s_beg,$s_end);
            $count ++;
        }
    }
    next if($count == 0); # hit no results.
    $maxCount=($maxCount>$count)?$maxCount:$count;
    my ($RNAL,$RNAR)=&MinMax(@RNAPos);
    my $gap5=$q_beg-$RNAL;
    my $gap3=$RNAR-$q_end;
    ($gap5,$gap3)=($gap3,$gap5) if($q_str eq "-");
    my $ids=join"\t",($q_id,@Rids);
    print ID $ids,"\n";
#   my $P1=join"\t",($q_id,"cov","len",$q_beg,$q_end,$q_str,@Rhits);
    my $P1=join"\t",($q_id,"cov","len",$FigL,$FigR,$q_str,$RNAL,$RNAR,$q_str,$q_beg,$q_end,$q_str,@Rhits);
    my $P2=$gap5."\t".$gap3;
    my $pos=join"\t",($q_beg,$q_end,$q_str,"");
    @{$Cakes{$q_id}}=($P1,$P2,$pos,$count);
}
close F;
close ID;

open OUT,"> $opt{o}" or die "Cannot open file $opt{o}:$!";
foreach my $p(sort keys %Cakes){
    my($P1,$P2,$pos,$count)=@{$Cakes{$p}};
    my $extra=$maxCount-$count;
    my $add=$pos x $extra;
    my $out=join"\t",($P1,$add,$P2);
    $out=~s/\t+/\t/g;
    print OUT $out,"\n";
}
close OUT;

my $columns=$maxCount * 3 + 14;
print "Max_hit_sRNAs & Columns_in_all:\t",$maxCount,"\t",$columns,"\n";

# Cal the overlap between two sequences
# , with the coordinations.
sub overlap{
    my($a_beg,$a_end, $b_beg,$b_end)=@_;
    my $a_len=$a_end-$a_beg+1;
    my $b_len=$b_end-$b_beg+1;
    my $t1=$a_end-$b_beg+1;
    my $t2=$b_end-$a_beg+1;
    my $over;
    if($a_beg<$b_end && $a_end>$b_beg){
        if($a_beg<$b_beg){
            $over=($t1>$b_len)?$b_len:$t1;
        }else{
            $over=($t2>$a_len)?$a_len:$t2;
        }
    }else{
        $over=0;
    }
    return $over;
}

sub MinMax{
    my ($min,$max)=@_;
    foreach my $i (@_){
        $min=($min<$i)?$min:$i;
        $max=($max>$i)?$max:$i;
    }
    return ($min,$max);
}
