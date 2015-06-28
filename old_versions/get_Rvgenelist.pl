#从*.gff文件提取Rv的基因列表。包含位置信息，正负链。

#! /user/perl -w
use strict;

open IN, $ARGV[0] or die "can't open the input file $!";
open INFO, $ARGV[1] or die "can't open the input file $!";
open OUT, ">",$ARGV[2] or die "can't open the input file $!";

my (%h,$cot);
while(<INFO>){      #读取Bacteria_info文件，以第4列 name为keys，第9列 为description.
	my @ps = split /\t/;
	$h{$ps[3]} = $ps[8];
$cot ++;	
	}

my ($num,$na);
while(<IN>){
	if($_ =~/RefSeq\sgene/){
		$_ =~ /locus\_tag\=(\w+);/;     #匹配上Rv0001；
		my $Rv = $1;
		if($_ =~ /ID\=NC_000962\.2\:/){
			$_ =~ /NC_000962\.2\:(\w+)\;/;
			$na = $1;
		}
		elsif($_ =~ /\;note\=tRNA/){
			$_ =~ /note\=(tRNA\-\w+)\%/;
			$na = $1;
			}
		else {$na = '-';}
				
		my @ps = split /\t/;
		print OUT "$Rv\t$na\t$ps[3]\t$ps[4]\t$ps[6]\t$h{$Rv}\n";
$num ++;				
	}else{next;}
			
}

print "Total gene\:$num\n";

close IN;
close INFO;
close OUT;