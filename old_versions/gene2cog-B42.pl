#########################################
#  Author: wangming  wangm08@yeah.net   #
#    					#
#  Version: 2010-07-15 	1.0		#
#					#
#########################################


#! /user/perl -w
use strict;
use Getopt::Long;

my ($help,$input1,$Head1,$input2,$Head2);
GetOptions(
	"help"=>\$help,	
	"clu=s"=>\$input1,
	"gene=s"=>\$input2,
	"Head1"=>\$Head1,
	"Head2"=>\$Head2,
);

if (!defined $input1 || !defined $input2 || defined $help) {
	&usage();
	exit 1;
}

open CLU, $input1 or die "Cannot open file $!.";
open GEN, $input2 or die "Cannot open file $!.";
open OUT, ">",$input2.'.cog.xls' or die "Cannot open file $!.";
open OUT2, ">",$input2.'.stat.xls' or die "Cannot open file $!.";

my (%clu,%gene,%unk,%des,$cont1,$cont2,$cont3);

#if (defined $Head1) {
#	<CLU>;
#	}

while(<CLU>){
	chomp;
	my @ps = split /\t/;
	$clu{$ps[1]} = {                	#读取第1列：序列，第3列： 数量，
			"A" => "$ps[0]",
			"Num" => "$ps[2]",		
		};
LINE:	my $tab = pop(@ps);		        #从尾部开始提取$ps[*]分割的单元，B42_refGL000021,COG0515 以这样一个单位，用tab分隔。
			if($tab =~/^B42\_refGL\d/){
				my @ges = split /,/, $tab;     	#分割 B42_refGL000080,COG2838 
				$clu{$ps[1]}->{$ges[0]} = 1; 	#设定两层hash，第一层keys为 功能类别，  第二层keys为 gene name。
				$unk{$ges[0]} = 1;
goto LINE;							#读取这一个function的所有基因。
				}
				else {next;}			#当一行中，Gene name提取完毕，转到读取下一个function类。
		}		

if (defined $Head2) {
	<GEN>;
	}
while(<GEN>){	#以gene列表的第一列gene name 为keys，将列表的一行存为values。
	chomp;
	my @ps = split /\t/;
	$gene{$ps[0]} = $_;	
$cont1 ++;	
	}
	
foreach my $c (sort keys %clu){					#以功能分类function 第一层hash ，遍历。。。
	print OUT $clu{$c}{'A'}."\t".$c."\n";			#print出每一个类别的信息，单独做一行。
	foreach my $g (sort keys %{$clu{$c}}){ 		# 遍历每一个功能分类下的所有基因。
		if(exists $gene{$g}){			#如果该基因在gene列表中存在，则打印出基因信息。
		print OUT $gene{$g}."\n";           
$cont2 ++;				
			}	else {next;}		#如果该基因在gene列表不存在，则转入下一个循环，判断下一个基因。
		}
print OUT2 "$clu{$c}{'A'}\t$c\t$cont2\n";		#在每一个function的所有基因遍历结束后，记录gene列表中在该function的数量，保存到 *.stat文件中。
$cont2 = 0;						#将$cont2归零，继续统计下一个function类的基因数量。
	}
	
print OUT "Not included\n";	
		
foreach my $ta (sort keys %gene)	{  		#最后打印出gene列表存在，但cluster中不存在的基因。命名为“Not included”.
	unless (exists $unk{$ta}) {
		print OUT "$gene{$ta}\n";
$cont3 ++;
		}
	}

print OUT2 "NO\tNot included\t$cont3\n";

close CLU;
close GEN;
close OUT;
close OUT2;


sub usage{
	print << "USAGE";
descripyion
Cluster file:
V       Defense mechanisms      42      B42_refGL000220,COG1132 B42_refGL000431,COG1680 ......
		
Gene file:
B42_refGL000053			read1	read2	……
B42_refGL000045			read1	read2 	……

example:	
	perl gene2cog-B42.pl  -h
	perl gene2cog-B42.pl  -clu B42.class.cog   -gene allgenelist.xls		 -Head1 0 	-Head2 1
output:
input.cog.xls	input.stat

USAGE
}
