#! /user/perl -w
use strict;
use Getopt::Long;

my ($help,$all_des,$all_exp,$input,$output,$Head);
GetOptions(
	"help"=>\$help,	
	"descrip=s"=>\$all_des,
	"express=s"=>\$all_exp, 
	"input=s"=>\$input,
	"output=s"=>\$output,
	"Head"=>\$Head,
);
if (!defined $input || defined $help || !defined $all_des || !defined $all_exp) {
	&usage();
	exit 1;
}

open DES, $all_des  or die "Cannot open file $!.";  #读取注释信息数据的总表。
open EXP, $all_exp  or die "Cannot open file $!."; #读取表达量数据的总表。

my (%hd, %tr,%he, $cont1, $cont2);
while(<DES>){
	chomp;
	my @ps = split /\t/;
	$hd{$ps[0]} = $_;
	$tr{$ps[6]} = $ps[0];
$cont1 ++;
	}
	
while(<EXP>){
	chomp;
	my @ps = split /\t/;
	$he{$ps[0]} = $_;
$cont2 ++;
	}

$output ||= $input;

open IN, $input or die "Cannot open file $!."; #读入需要注释的基因列表。
open OUTd, ">",$output.'-des.xls' or die "Cannot open file $!.";
open OUTe, ">",$output.'-exp.xls' or die "Cannot open file $!.";
if (defined $Head) {
	<IN>;
}

my $cont3;
while(<IN>){
	chomp;
	my @ps = split /\t/;
	
	if($ps[0] =~/^B42/){
	print OUTd "$hd{$ps[0]}\n";
	print OUTe "$he{$ps[0]}\n";
$cont3 ++;
	  }
	  elsif($ps[0] =~/^Rv/){
	  	   	print OUTd "$hd{$tr{$ps[0]}}\n";
	       print OUTe "$he{$tr{$ps[0]}}\n";
$cont3 ++;
	  	}
}	
print "DES\:$cont1\nEXP\:$cont2\nIN\:$cont3\n";

close DES;
close EXP;
close IN;
close OUTd;
close OUTe;

sub usage{
	print << "USAGE";
name
	get_des.pl
descripyion
	get gene description and expression according to the geneID.
	-input  (str) the input file \$1 B42_refGL000001 or Rv0001 ;
	-output (str) output 
	-Head        the input with head ( 1 line)
	-help         help
author

version
	1.0    2010-6-28 10:40
example	
	perl get_des.pl -h
	perl get_des.pl -descrip all_des.xls -express all_exp.xls  -input gene-list.xls  -Head  1
USAGE
}
