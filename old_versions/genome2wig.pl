#!/usr/bin/perl
$f=shift;
$strand=shift;
$name=shift;
$des=shift;
#print "$strand\n";
unless(-e $f){print "Input match_genome.txt file!\nperl genome2wig.pl  match_genome.txt  +/-  name description\n";exit;}
#if($strand ne "p" or $strand ne "m"){print "Strand must be +/-!\nperl bam2wig.pl match_genome.txt +/- name description\n";exit;}

if($strand eq "p"){$strand="+";}
if($strand eq "m"){$strand="-";}

print "track type=wiggle_0 name=$name description=$des\n";
open IN,$f;
while(<IN>){
	chomp;
	@d=split(/\t/,$_);
	next if($d[4] ne $strand);
	
	$n=$d[6];	#Read count.
	foreach $i($d[2]..$d[3]){
		$count{$d[1]}{$i}+=$n;
		$chrome{$d[1]}=1;
		}
	}
close IN;

foreach $chr_id(sort keys %chrome){
	print "variableStep chrom=$chr_id\n";
	foreach $loc(sort { $a<=>$b } keys %{$count{$chr_id}}){
		print "$loc\t$count{$chr_id}{$loc}\n";
		}
	}
exit;