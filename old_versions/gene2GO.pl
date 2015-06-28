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
open CLU, $input1 or die "Cannot open file $input1";        #读取cluster文件。
open GEN, $input2 or die "Cannot open file $input2";         #读取基因信息，或是 基因列表，包含完整的基因表达or注释信息。
open Rv, '/share/raid12/wangming/reference_genomes/M.tb/Rv_genelist.xls' or die "Cannot open file $!.";         #读取H37Rv中基因起始位置, 作为keys 保存为$Rv{$begin} = $name;
open OUT, ">",$input1.'.cluster.xls' or die "Cannot open file $!.";    #输出文件1，根据cluster中，每一个类中基因的详细信息。
open OUT2, ">",$input1.'.stat.xls' or die "Cannot open file $!.";     #统计cluster中，基因总数。

#=================================
#   按照cluster的格式保存为 %clu.
#=================================

my (%clu,%h,%gene,%con,%all,%Rv,$cont1,$cont2,$cont3);

while(<Rv>){
	chomp;
	my @ps = split /\t/;
	$Rv{$ps[2]} = $ps[0]; 
}

while(<CLU>){
unless($_ =~/B42\w+|Rv\w+/){
	print OUT $_;     #写出给每一个cluster的类的标题
	print OUT2 $_;
	next;
	}
	chomp;
	my @ps = split /\t/;
	my $list = pop(@ps);
	my $line = join "\t", @ps;     #将每一行除基因之外的部分，单独写出。
		$h{$ps[0]} = $line;
		my @genes = split /\,/, $list; 		#分割最后一列，B42_refGL001434,B42_refGL001436,B42_refGL001433
LINE:		my $g = shift (@genes);		#从数组开头提取基因名称 Rv0001。。。
				unless(defined $g){ next;}
  if($g =~ /NC\_000962/){
	   if($g =~ /\:c\d+\-(\d+)/) {$g= $Rv{$1};} #$1等于该基因起始位点，转换$nt 为Rv0000.
		if($g =~ /\:(\d+)\-\d+/) {$g = $Rv{$1};} 
	  }	
				$clu{$ps[0]}->{$g} = 1;	#给存在的基因名称赋值为1.
				$con{$ps[0]} ++;		          #统计每一个function 类中基因数目。
				$all{$g} ++;						#记录cluster类中每一个基因出现的次数。
		goto LINE;		
}		

while(<GEN>){
	unless($_ =~/^B42\_refGL|^Rv\d+/){ next;}  #匹配第一列是否为基因名称，B42_refGL*  or  Rv0000
	chomp;
	my @ps = split /\t/;
	my $g = shift(@ps);
	my $line = join "\t", @ps;
	$gene{$g} = $line;				#将gene name 作为keys，第2列之后的信息为 values  存入%gene。
$cont1 ++;	
	}
	
#=================================	
#	根据cluster的列表格式，给每一类的基因，添加信息。	
#=================================	

foreach my $c (keys %clu){					#以功能分类function 进行循环。。。
	print OUT "$h{$c}\t$con{$c}\n";			#print出每一个类别的信息，单独做一行，尾部附上这一类别的基因数。
	foreach my $g (keys %{$clu{$c}}){ 		# 对每一个功能分类下的所有基因 进行循环。
		unless(exists $gene{$g}){next;}			#如果该基因在gene列表中存在，则打印出基因信息，否则判断下一个基因。
		print OUT "$g\t$all{$g}\t$gene{$g}\n";   #将该基因的INFO，及在cluster中出现的次数加入。
$cont2 ++;				
		  }
print OUT2 "$h{$c}\t$cont2\n";		#记录每一个function类的信息及基因数目，到 *.stat文件中。
$cont2 = 0;						#将$cont2归零，继续统计下一个function类的基因数量。
	}
	
print OUT "Not included\n";	
foreach my $g (keys %gene)	{  	#最后打印出gene列表存在，但cluster中不存在的基因。命名为“Not included”.
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
