#=======================================
#  Author: wangming  wangm08@yeah.net                           
#    					                                                                 
#  Version: 2010-07-17 	1.0		                                          
#					                                                                    
#=======================================

#! /user/perl -w
use strict;

unless(defined @ARGV) {
	&usage();
	exit 1;
}
my $input1 = $ARGV[0];
my $input2 = $ARGV[1];
open CLU, $input1 or die "Cannot open file $input1";        #��ȡcluster�ļ���
open GEN, $input2 or die "Cannot open file $input2";         #��ȡ������Ϣ������ �����б����������Ļ�����orע����Ϣ��
open Rv, '/share/raid12/wangming/reference_genomes/M.tb/Rv_genelist.xls' or die "Cannot open file $!.";         #��ȡH37Rv�л�����ʼλ��, ��Ϊkeys ����Ϊ$Rv{$begin} = $name;
open OUT, ">",$input1.'.cluster.xls' or die "Cannot open file $!.";    #����ļ�1������cluster�У�ÿһ�����л������ϸ��Ϣ��
open OUT2, ">",$input1.'.stat.xls' or die "Cannot open file $!.";     #ͳ��cluster�У�����������

#=================================
#   ����cluster�ĸ�ʽ����Ϊ %clu.
#=================================

my (%clu,%h,%gene,%con,%all,%Rv,$cont1,$cont2,$cont3);

while(<Rv>){
	chomp;
	my @ps = split /\t/;
	$Rv{$ps[2]} = $ps[0]; 
}

while(<CLU>){
unless($_ =~/B42\w+|Rv\w+/){
	print OUT $_;     #д����ÿһ��cluster����ı���
	print OUT2 $_;
	next;
	}
	chomp;
	my @ps = split /\t/;
	my $list = pop(@ps);
	my $line = join "\t", @ps;     #��ÿһ�г�����֮��Ĳ��֣�����д����
		$h{$ps[0]} = $line;
		my @genes = split /\,/, $list; 		#�ָ����һ�У�B42_refGL001434,B42_refGL001436,B42_refGL001433
LINE:		my $g = shift (@genes);		#�����鿪ͷ��ȡ�������� Rv0001������
				unless(defined $g){ next;}
  if($g =~ /NC\_000962/){
	   if($g =~ /\:c\d+\-(\d+)/) {$g= $Rv{$1};} #$1���ڸû�����ʼλ�㣬ת��$nt ΪRv0000.
		if($g =~ /\:(\d+)\-\d+/) {$g = $Rv{$1};} 
	  }	
				$clu{$ps[0]}->{$g} = 1;	#�����ڵĻ������Ƹ�ֵΪ1.
				$con{$ps[0]} ++;		          #ͳ��ÿһ��function ���л�����Ŀ��
				$all{$g} ++;						#��¼cluster����ÿһ��������ֵĴ�����
		goto LINE;		
}		

while(<GEN>){
	unless($_ =~/^B42\_refGL|^Rv\d+/){ next;}  #ƥ���һ���Ƿ�Ϊ�������ƣ�B42_refGL*  or  Rv0000
	chomp;
	my @ps = split /\t/;
	my $g = shift(@ps);
	my $line = join "\t", @ps;
	$gene{$g} = $line;				#��gene name ��Ϊkeys����2��֮�����ϢΪ values  ����%gene��
$cont1 ++;	
	}
	
#=================================	
#	����cluster���б��ʽ����ÿһ��Ļ��������Ϣ��	
#=================================	

foreach my $c (keys %clu){					#�Թ��ܷ���function ����ѭ��������
	print OUT "$h{$c}\t$con{$c}\n";			#print��ÿһ��������Ϣ��������һ�У�β��������һ���Ļ�������
	foreach my $g (keys %{$clu{$c}}){ 		# ��ÿһ�����ܷ����µ����л��� ����ѭ����
		unless(exists $gene{$g}){next;}			#����û�����gene�б��д��ڣ����ӡ��������Ϣ�������ж���һ������
		print OUT "$g\t$all{$g}\t$gene{$g}\n";   #���û����INFO������cluster�г��ֵĴ������롣
$cont2 ++;				
		  }
print OUT2 "$h{$c}\t$cont2\n";		#��¼ÿһ��function�����Ϣ��������Ŀ���� *.stat�ļ��С�
$cont2 = 0;						#��$cont2���㣬����ͳ����һ��function��Ļ���������
	}
	
print OUT "Not included\n";	
foreach my $g (keys %gene)	{  	#����ӡ��gene�б���ڣ���cluster�в����ڵĻ�������Ϊ��Not included��.
	unless (exists $all{$g}) {
		print OUT "$g\t0\t$gene{$g}\n";
$cont3 ++;
		}
	}
print OUT2 "NO\tNot included\t$cont3\n";


sub usage{
	print << "USAGE";
perl gene2clu.pl  cluster.txt   allgenelist.xls	
output:
input.cog.xls	input.stat
USAGE
}
