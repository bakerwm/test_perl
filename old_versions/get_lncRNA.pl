### !-------------------------------------------------------
### delete, replaced by: search_cov_tags.pl


#!/usr/bin/perl -w
use strict;
use warnings;

my $infile=shift;
my $reffile=shift;
my $strain=shift;
my $fd=shift;
unless(-e $infile or -e $reffile or -e $strain or -e $fd){
	print "perl $0 PosCoverage.txt reference_file strain fold_change\n";
	exit;
	}
my %Pos_Cov;
open IN,"$infile" || die "$!\n";
while(<IN>){
	chomp;
	my @d=split;
	$Pos_Cov{$d[1]}=$d[2];
}
close IN;

my $len_ref;
open IN,"$reffile";
while(<IN>){
	chomp;
	next if(/^>/);
	$len_ref+=length $_;
}

print "ID\tchain\tBegin:Coverage\tMax:coverage\tEnd:coverage\tLength\n";
my $num=0;
my $tem=1;
my $lase_seed_end=0;
for(my $i=1;$i<$len_ref;$i=$tem){
	$tem=$i+1;
### find the highest position ###
	next if(not defined $Pos_Cov{$i});
	next if($Pos_Cov{$i} < 100);
	
	my $maxpos=$i;
	my $maxcov=$Pos_Cov{$i};
	my $start_pos=$i;
	for(my $j=$i;$j<$i+1000;$j++){
		last unless($Pos_Cov{$j} or $Pos_Cov{$j+1});
		if($Pos_Cov{$j} >= $maxcov){
			$maxcov=$Pos_Cov{$j};
			$maxpos=$j;
		next unless($Pos_Cov{$j+1});
			if($Pos_Cov{$j+1} < $maxcov){
				my $flag=0;
				for my $l(1..20){
					next if(not defined $Pos_Cov{$j+$l});
					if(($Pos_Cov{$j+$l} > $Pos_Cov{$j})){
						$flag=1;
					}
				}
			 	last if($flag == 0);
			}
		}
	}
	

	my $halfmaxcov=int($maxcov*$fd);  ### define the limit expression level
	my $leftpos=$maxpos;
	my $rightpos=$maxpos;
	my $l_s;
	
	if($lase_seed_end > $maxpos-1000){$l_s=$lase_seed_end+1;}else{$l_s=$maxpos-1000;}
	for(my $k=$maxpos;$k >= $l_s;$k--){
#		if($k=$lase_seed_end+1){$leftpos=$k;last;}
		unless($Pos_Cov{$k-1}){$leftpos=$k;last;}
		if(($Pos_Cov{$k} >= $halfmaxcov) && ($Pos_Cov{$k-1} < $halfmaxcov)){
			$leftpos=$k;	
			last;
		}
	}
	
	for(my $m=$maxpos;$m <= $maxpos+1000;++$m){
		unless($Pos_Cov{$m+1}){$rightpos=$m;last;}
		if(($Pos_Cov{$m} >= $halfmaxcov) && ($Pos_Cov{$m+1} < $halfmaxcov)){
			$rightpos=$m;
			last;
		}
	}
=pod
	my 	$len=$rightpos-$leftpos+1;
	$tem=$rightpos+1;
	$lase_seed_end=$rightpos;
	next if($len<16);
	++$num;
	$num=sprintf "%.4d",$num;
	print "NUM$num\t$strain\t$leftpos:$Pos_Cov{$leftpos}\t$maxpos:$maxcov\t$rightpos:$Pos_Cov{$rightpos}\t$len\n";
=cut		
#	next if($leftpos-$start_pos<50);
	my $temp_max_pos=0;my $temp_max_cov=0;
	
	for(my $temp=$start_pos;$temp<$leftpos;$temp++){
		next unless($Pos_Cov{$temp});
		if($Pos_Cov{$temp}>$temp_max_cov){$temp_max_pos=$temp;$temp_max_cov=$Pos_Cov{$temp};}}

	my $left_temp_pos=$temp_max_pos;my $right_temp_pos=$temp_max_pos;
	$halfmaxcov=int($temp_max_cov*$fd);
	for(my $s=$temp_max_pos;$s>=$start_pos;$s--){
		unless($Pos_Cov{$s}){$left_temp_pos=$s+1;last;}
		if($Pos_Cov{$s} >= $halfmaxcov){$left_temp_pos=$s;}
		else{$left_temp_pos=$s-1;last;}}

	for(my $r=$temp_max_pos;$r<$leftpos;$r++){
		unless($Pos_Cov{$r}){$right_temp_pos=$r-1;last;}
		if($Pos_Cov{$r} >= $halfmaxcov ){$right_temp_pos=$r;}
		else{$right_temp_pos=$r-1;last;}	}

	my 	$len_main=$rightpos-$leftpos+1;
	$tem=$rightpos+1;
	$lase_seed_end=$rightpos;
	if($len_main<16){
		my $len=$right_temp_pos-$left_temp_pos+1;
		next if($len<16);
		++$num;
		$num=sprintf "%.4d",$num;
		print "NUM$num\t$strain\t$left_temp_pos:$Pos_Cov{$left_temp_pos}\t$temp_max_pos:$temp_max_cov\t$right_temp_pos:$Pos_Cov{$right_temp_pos}\t$len\n";	
		next;
		}
	++$num;
	$num=sprintf "%.4d",$num;
	print "NUM$num\t$strain\t$leftpos:$Pos_Cov{$leftpos}\t$maxpos:$maxcov\t$rightpos:$Pos_Cov{$rightpos}\t$len_main\n";
}
exit;

