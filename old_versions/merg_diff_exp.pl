#�ϲ������������б�����hash�洢�� $h{$�б�}{$gene}=[split], ����Ĺ�ϣ��
#���жϻ����б��У�ĳЩ�����Ƿ���ڣ�������p-value��FDR�������ڣ����á�-�����档
#

#��/user/bin/perl -w
use strict;

unless (defined @ARGV) {
	&usage();
	exit 1;
}

my $des = shift @ARGV;   #��ȡ����ע���ܱ� all_des.xls
open Rv, $des or die "Cannot open file  $!.";  

my (%Rv_pos,%all,%all_gene);
my ($cont1,$num1,$num2,$num3);
while(<Rv>){
	chomp;
	my @ps = split /\t/;
	$all_gene{$ps[0]} = [@ps]; #�������еĻ������ƣ���1�У�B42_refGL000000������Rv0000��
	next unless /^(\Rv\d+)\s/; #ƥ����Rv��ͷ���С�
	$Rv_pos{$ps[9]} = $ps[6];     #��beginΪkey��Rv0000Ϊֵ���� $hr��
$cont1 ++;
	}

my $diff_1 = shift @ARGV;
my $diff_2 = shift @ARGV;
my $diff_3 = shift @ARGV;	
my $expr   = shift @ARGV;

open Diff_1, $diff_1 or die "Cannot open file $!.";  #Diff_1 ��ָ�������б� Rv-0vs42-0-pvalue.xls
open Diff_2, $diff_2 or die "Cannot open file $!.";  #Diff_2 ��ָ�������б� 42-0vs42-6-pvalue.xls
open Diff_3, $diff_3 or die "Cannot open file $!.";  #Diff_3  ��ָ�������б� Rv-0vs42-6-pvalue.xls	
open OUT, '>',$expr or die "Cannot open file $!."; #���ɱ������Ϣ�ܱ�



while(<Diff_1>){
	next if $_ =~ /^geneID/; #ȥ����ͷ����geneIDΪ��ͷ��һ�С�
	chomp;
	my @ps = split /\t/;
	my $g = shift @ps;
	if($g =~ /NC\_000962/){
		my $begin;
		if($g =~ /\|\:\w+\-(\d+)/){$begin = $1;} #��Rv��ʼλ�����Ϣת��Ϊ������Rv0000.
		if($g =~ /\|\:(\d+)\-/){$begin = $1;}
		$g = $Rv_pos{$begin};	 #����ʼλ�����Ϣת��Ϊ�������ơ�
	 } 
	$all{Diff_1}{$g} = [@ps]; #������ȥ��ԭ��1��-geneID����ʣ�µ�2��---����С�
$num1 ++;		
}

while(<Diff_2>){
	next if $_ =~ /^geneID/; #ȥ����ͷ����geneIDΪ��ͷ��һ�С�
	chomp;
	my @ps = split /\t/;
	my $g = shift @ps;
	if($g =~ /NC\_000962/){
		my $begin;
		if($g =~ /\|\:\w+\-(\d+)/){$begin = $1;} #��Rv��ʼλ�����Ϣת��Ϊ������Rv0000.
		if($g =~ /\|\:(\d+)\-/){$begin = $1;}
		$g = $Rv_pos{$begin};	 #����ʼλ�����Ϣת��Ϊ�������ơ�
	 } 
	$all{Diff_2}{$g} = [@ps]; #������ȥ��ԭ��1��-geneID����ʣ�µ�2��---����С�
$num2 ++;		
}

while(<Diff_3>){
	next if $_ =~ /^geneID/; #ȥ����ͷ����geneIDΪ��ͷ��һ�С�
	chomp;
	my @ps = split /\t/;
	my $g = shift @ps;
	if($g =~ /NC\_000962/){
		my $begin;
		if($g =~ /\|\:\w+\-(\d+)/){$begin = $1;} #��Rv��ʼλ�����Ϣת��Ϊ������Rv0000.
		if($g =~ /\|\:(\d+)\-/){$begin = $1;}
		$g = $Rv_pos{$begin};	 #����ʼλ�����Ϣת��Ϊ�������ơ�
	 } 
	$all{Diff_3}{$g} = [@ps]; #������ȥ��ԭ��1��-geneID����ʣ�µ�2��---����С�
$num3 ++;		
}

foreach my $g (sort keys %all_gene){
#	my $null = "\t\-" x 9;

if(!exists $all{Diff_1}{$g}){
	$all{Diff_1}{$g} = [qw(0 0 0 0 0 - - - - )];	
	}	
if(!exists $all{Diff_2}{$g}){
	$all{Diff_2}{$g} = [qw(0 0 0 0 0 - - - - )];
	}
if(!exists $all{Diff_3}{$g}){
	$all{Diff_3}{$g} = [qw(0 0 0 0 0 - - - - )];
	}	

my $Read = $all{Diff_1}{$g}[1]."\t".$all{Diff_2}{$g}[1]."\t".$all{Diff_3}{$g}[2];
my $RPKM = $all{Diff_1}{$g}[3]."\t".$all{Diff_2}{$g}[3]."\t".$all{Diff_3}{$g}[4];
my $log2  = $all{Diff_1}{$g}[5]."\t".$all{Diff_2}{$g}[5]."\t".$all{Diff_3}{$g}[5];
#my $p_val = $all{Diff_1}{$g}[7]."\t".$all{Diff_2}{$g}[7]."\t".$all{Diff_3}{$g}[7];
my $FDR  = $all{Diff_1}{$g}[8]."\t".$all{Diff_2}{$g}[8]."\t".$all{Diff_3}{$g}[8];	
	
print OUT $g."\t".$Read."\t".$RPKM."\t".$log2."\t".$FDR."\n";
	
}

close Rv;
close Diff_1;
close Diff_2;
close Diff_3;
close OUT;


sub usage{
	print << "USAGE";
perl merg_diff_exp.pl  all_des.xls  diff_exp1.xls  diff_exp2.xls  diff_exp3.xls  all_exp.xls
USAGE
}
