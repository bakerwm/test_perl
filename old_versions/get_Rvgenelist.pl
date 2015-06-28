#��*.gff�ļ���ȡRv�Ļ����б�����λ����Ϣ����������

#! /user/perl -w
use strict;

open IN, $ARGV[0] or die "can't open the input file $!";
open INFO, $ARGV[1] or die "can't open the input file $!";
open OUT, ">",$ARGV[2] or die "can't open the input file $!";

my (%h,$cot);
while(<INFO>){      #��ȡBacteria_info�ļ����Ե�4�� nameΪkeys����9�� Ϊdescription.
	my @ps = split /\t/;
	$h{$ps[3]} = $ps[8];
$cot ++;	
	}

my ($num,$na);
while(<IN>){
	if($_ =~/RefSeq\sgene/){
		$_ =~ /locus\_tag\=(\w+);/;     #ƥ����Rv0001��
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