#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;

sub help{
        print STDERR <<EOF;
Usage: perl  $0  -i coverage  -s strand  -c cutoff  -o outfile
Options:
    -i <file>       : coverage file, 3-line tab-delimited file in the following format:
                      <Genome_ID> <Position> <Coverage>
                      Position should be continuous from 1 to end.
    -s [+/-]        : The strand.
    -c [0-1]        : The cutoff [0-1] to determin the edge. Default 0.5 .
    -o <file>       : The result file.
EOF
}

my ($infile, $strand, $cutoff, $outfile, $help);
&GetOptions(
            "i=s"   =>  \$infile,
            "s=s"   =>  \$strand,
            "c=s"   =>  \$cutoff,
            "o=s"   =>  \$outfile,
            "h"     =>  \$help,
            );
$cutoff = 0.5 if(not $cutoff);
if($help or not $infile or not $strand or not $outfile){
    &help();
    exit(1);
}

my $label = ($strand eq '+')?'_p':'_n';

#####################
# Read coverage file
#####################
my %h_cov;
my $genome_end = 0;
my %mark;
open F,"< $infile" or die;
my $check_pos = 0;
while(<F>){
    chomp;
    my @ps = split(/\t/);

if(($ps[1] - $check_pos) > 1){
print "Warn: check the input coverage file: $infile
The positions in line 2 are not continuous from 1 to end.\n";
exit(1);    
    }
    $check_pos = $ps[1];

    $h_cov{$ps[1]} = $ps[2];    #
    $genome_end    = $ps[1];    # Length of genome.
    $mark{$ps[1]}  = 0;         # mark the sRNAs.
}
close F;

# Print Header.
open OUT,">$outfile" or die;
my $Header = join"\t",("ID", "Strand", "Begin\:Cov", "Max\:Cov", "End\:Cov", "Length");
print OUT $Header,"\n";

#######################
# Searching sRNAs 
#######################
#my ($start, $end)         = ("1", $genome_end);
my ($gap_start, $gap_end) = ("1", $genome_end);
my $id_count = 1;

for(my $n=1; $n <= $genome_end; $n++){
    if($mark{$n} == 1){
        next;
    }else{
# Determin the range of gap.
FIN:
        $gap_start = $n;
        for(my $p=$gap_start; $p<=$genome_end; $p++){
            if($mark{$p} == 1){
                $gap_end = $p - 1;
                last;
            }else{
                next;
            }
        }
        $gap_end = $genome_end if($gap_start > $gap_end);
        
        my ($begin_pos, $begin_cov, $max_pos, $max_cov, $end_pos, $end_cov) = &get_sRNA($gap_start, $gap_end);
        my $id     = sprintf "%0.4d",$id_count;
        my $length = $end_pos - $begin_pos + 1;
        my $out    = join"\t",("NUM$id$label",$strand,"$begin_pos\:$begin_cov","$max_pos\:$max_cov", "$end_pos\:$end_cov", $length);

        for($begin_pos..$end_pos){
            $mark{$_} = 1;
        }
        
        if($max_cov >= 100 && $length >= 16){
            print OUT $out,"\n";
            $id_count ++;
        }
        
        if($begin_pos > $gap_start){
            goto FIN;
        }else{
            next;
        }
    }
}

########################################
# Subroutine:
# Find one sRNAs from a specific region;
########################################
sub get_sRNA{
    my ($start_pos, $end_pos) = @_;
    my ($start_cov, $end_cov) = ($h_cov{$start_pos}, $h_cov{$end_pos});
    my ($max_pos,   $max_cov) = ($start_pos, $start_cov);
# Scanning for the max position.
    for(my $i=$start_pos; $i<=$end_pos; $i++){
        my $check_cov = $h_cov{$i};
        if($check_cov >= $max_cov){
            ($max_pos, $max_cov) = ($i, $check_cov);
        }    
    }
# Scanning for the 3' end.
    for(my $j=$max_pos; $j<=$end_pos; $j++){
        my $check_cov = $h_cov{$j};
        if($check_cov < ($max_cov * $cutoff)){
            $end_pos = $j - 1;
            $end_cov = $h_cov{$end_pos};
            last;
        }else{
            next;
        }
    }
# search toward the 5' end.
    for(my $k=$max_pos; $k>=$start_pos; $k--){
        my $check_cov_5end = $h_cov{$k};
        if($check_cov_5end < ($max_cov * $cutoff)){
            $start_pos = $k + 1;
            $start_cov = $h_cov{$start_pos};
            last;
        }else{
            next;
        }
    }
# return 
    return($start_pos, $start_cov, $max_pos, $max_cov, $end_pos, $end_cov);
}
