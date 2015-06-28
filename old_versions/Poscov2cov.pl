#!/usr/bin/perl  -w

use strict;
use Data::Dumper;
use POSIX qw(ceil floor);

my $usage = 'Usage: perl  Poscov2cov.pl  Poscovage  input_list  > output_file
##Format for Poscovage
1           2
position    cov

##Format for input_list
id  exp len begin   end strand
--  --  --  --      --  --';

my $Poscov     = shift or die $usage;
my $input_list = shift or die $usage;

my (
    %ha,
    %stat,
    );

open F,$Poscov or die;
while(<F>){
    chomp;
    my @ps=split/\t/;
    $ha{$ps[2]}->{$ps[0]} = $ps[1]; # key1: strand, key2: position.
}
close F;

my $gene_number = 0;
open F,$input_list or die;
while(<F>){
$gene_number ++;    
    chomp;
    my @ps=split/\t/;
    my ($temp1, $temp2, $mR_len, $mR_begin, $mR_end, $mR_strand) = @ps;

    for(my $i = $mR_begin; $i <= $mR_end; $i ++){
        my $stat_pos;
        my $km;
        if($mR_strand eq "+"){        
           $km = sprintf "%.2f", (($i - $mR_begin + 1)/$mR_len * 100);
        }else{
           $km = sprintf "%.2f", (($mR_end - $i + 1)/$mR_len * 100);
        }
        $stat_pos = ceil($km);
        my $exp = (exists $ha{$mR_strand}->{$i})?$ha{$mR_strand}->{$i}:0;

        $stat{$stat_pos} += $exp;
    }
}
close F;

foreach(sort{$a<=>$b} keys %stat){
    my $exp = sprintf "%.2f", $stat{$_}/$gene_number;
    print $_,"\t",$exp,"\n";
}
