#!/usr/bin/perl -w
use strict;
use warnings;


my $usage="perl  PttRnt2table.pl   NC_000962.ptt  NC_000962.rnt  H37Rv";
my $PTT=shift or die"Input the *.ptt file: \n$usage";
my $rnt=shift or die"Input the *.rnt file: \n$usage";
my $Name=shift or die"Input the name of strain: \n$usage";

open F,$PTT or die;
my $F_out="";
while(<F>){
    next unless(/^\d+\.\./);
    chomp;
    my @p=split/\t/;
    my ($beg,$end)=($p[0]=~/(\d+)\.\.(\d+)/);
    my $len=$end-$beg+1;
    my $out=join"\t",($p[5],$p[4],$len,$beg,$end,$p[1],"mRNA",$p[8],"\n");
    $F_out.=$out;
}
close F;

open F,$rnt or die;
while(<F>){
    next unless(/^\d+\.\./);
    chomp;
    my @t=split/\t/;
    my ($beg,$end)=($t[0]=~/(\d+)\.\.(\d+)/);
    my $len=$end-$beg+1;
    my $tag=(/ribosomal RNA/)?"rRNA":((/Anticodon/)?"tRNA":"other");
    next if($tag eq "other");
    my $out=join"\t",($t[5],$t[4],$len,$beg,$end,$t[1],$tag,$t[8],"\n");
    $F_out.=$out;
}
close F;

my $outfile=$Name."_mtrRNA.table";
open OUT,">$outfile" or die;
print OUT $F_out;
close OUT;

