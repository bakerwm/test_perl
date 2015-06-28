#从B42.faa.H37RV.gbk.pep.reciprocal.best.list获取共有的基因列表，然后分别读取Rv和B42的基因列表。
#获得基因列表的总表。

#! /user/perl -w
use strict;

if (!defined @ARGV) {
	&usage();
	exit 1;
}

open BEST, $ARGV[0] or die "can't open the input file $!";
open B42, $ARGV[1] or die "can't open the input file $!";
open Rv, $ARGV[2] or die "can't open the input file $!";
open OUT, ">",$ARGV[3] or die "can't open the input file $!";

my ($num1,$num2,$num3,$num4);
my (%ha,%hb,%h2,%h3);

while(<BEST>){
	chomp;
	my @ps = split /\t/;
	$ha{$ps[0]} = {
#			'Rv' => $ps[5],
			'p-val' => $ps[10],
			'Iden' =>  $ps[11],
		};
	$hb{$ps[5]} = {
			'B42' => $ps[0],
			'p-val' => $ps[10],
			'Iden' =>  $ps[11],		
		};
	}
	
while(<B42>)	{
	chomp;
	my @ps = split /\t/;
	my $len = $ps[2] - $ps[1] + 1;
	$h2{$ps[0]} = {
		'len' => $len,
		'begin' => $ps[1],
		'end' => $ps[2],
		'str' => $ps[3],		
		};
}

while(<Rv>)	{
	chomp;
	my @ps = split /\t/;
	my $len = $ps[3] - $ps[2] + 1;
	$h3{$ps[0]} = {
		'len' => $len,
		'name' => $ps[1],
		'begin' => $ps[2],
		'end' => $ps[3],
		'str' => $ps[4],
		'des'	 => $ps[5],
		};
$num4 ++;		
}

	
foreach my $g (sort keys %h3){	
   if(exists $hb{$g}){
my $r = $hb{$g}{'B42'};  	
        my $part1 = "$r\tNA\t$h2{$r}{'len'}\t$h2{$r}{'begin'}\t$h2{$r}{'end'}\t$h2{$r}{'str'}\t";	
   	     my $part2 = "$g\t$h3{$g}{'name'}\t$h3{$g}{'len'}\t$h3{$g}{'begin'}\t$h3{$g}{'end'}\t$h3{$g}{'str'}\t";
  	     my $part3 = "$hb{$g}{'p-val'}\t$hb{$g}{'Iden'}\t";   	
   	print OUT $part1.$part2.$part3."$h3{$g}{'des'}\n"; 
$num1 ++;   	
      } else {
   	       my $part1 = "$g\tRS\t".'-	-	-	-	';
   	       my $part2 = "$g\t$h3{$g}{'name'}\t$h3{$g}{'len'}\t$h3{$g}{'begin'}\t$h3{$g}{'end'}\t$h3{$g}{'str'}\t";
    	   my  $part3 = '-	-	';
   	     print OUT $part1.$part2.$part3."$h3{$g}{'des'}\n";  
$num3 ++;   		
   	}	
}

foreach my $g (sort keys %h2){	
   if(exists $ha{$g}){
next; 	
  } else {   	
   my $part1 = "$g\tBS\t$h2{$g}{'len'}\t$h2{$g}{'begin'}\t$h2{$g}{'end'}\t$h2{$g}{'str'}\t";
   	my $part2 = "\-\t" x 6;
   	my $part3 = "\-\t" x 2;
   	print OUT $part1.$part2.$part3."\-\n";  
$num2 ++;   		
   	}	
}

print "both\:$num1\nB42\-spec\:$num2\nRv\-spec\:$num3\ntest\:$num4\n";

close BEST;
close B42;
close Rv;
close OUT;

sub usage{
	print << "USAGE";
perl specific.pl  best.xls 	B42_genelist.xls     Rv_genelist.xls    out.xls
USAGE
}