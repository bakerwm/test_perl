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
	$clu{$ps[1]} = {                	#��ȡ��1�У����У���3�У� ������
			"A" => "$ps[0]",
			"Num" => "$ps[2]",		
		};
LINE:	my $tab = pop(@ps);		        #��β����ʼ��ȡ$ps[*]�ָ�ĵ�Ԫ��B42_refGL000021,COG0515 ������һ����λ����tab�ָ���
			if($tab =~/^B42\_refGL\d/){
				my @ges = split /,/, $tab;     	#�ָ� B42_refGL000080,COG2838 
				$clu{$ps[1]}->{$ges[0]} = 1; 	#�趨����hash����һ��keysΪ �������  �ڶ���keysΪ gene name��
				$unk{$ges[0]} = 1;
goto LINE;							#��ȡ��һ��function�����л���
				}
				else {next;}			#��һ���У�Gene name��ȡ��ϣ�ת����ȡ��һ��function�ࡣ
		}		

if (defined $Head2) {
	<GEN>;
	}
while(<GEN>){	#��gene�б�ĵ�һ��gene name Ϊkeys�����б��һ�д�Ϊvalues��
	chomp;
	my @ps = split /\t/;
	$gene{$ps[0]} = $_;	
$cont1 ++;	
	}
	
foreach my $c (sort keys %clu){					#�Թ��ܷ���function ��һ��hash ������������
	print OUT $clu{$c}{'A'}."\t".$c."\n";			#print��ÿһ��������Ϣ��������һ�С�
	foreach my $g (sort keys %{$clu{$c}}){ 		# ����ÿһ�����ܷ����µ����л���
		if(exists $gene{$g}){			#����û�����gene�б��д��ڣ����ӡ��������Ϣ��
		print OUT $gene{$g}."\n";           
$cont2 ++;				
			}	else {next;}		#����û�����gene�б����ڣ���ת����һ��ѭ�����ж���һ������
		}
print OUT2 "$clu{$c}{'A'}\t$c\t$cont2\n";		#��ÿһ��function�����л�����������󣬼�¼gene�б����ڸ�function�����������浽 *.stat�ļ��С�
$cont2 = 0;						#��$cont2���㣬����ͳ����һ��function��Ļ���������
	}
	
print OUT "Not included\n";	
		
foreach my $ta (sort keys %gene)	{  		#����ӡ��gene�б���ڣ���cluster�в����ڵĻ�������Ϊ��Not included��.
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
B42_refGL000053			read1	read2	����
B42_refGL000045			read1	read2 	����

example:	
	perl gene2cog-B42.pl  -h
	perl gene2cog-B42.pl  -clu B42.class.cog   -gene allgenelist.xls		 -Head1 0 	-Head2 1
output:
input.cog.xls	input.stat

USAGE
}
