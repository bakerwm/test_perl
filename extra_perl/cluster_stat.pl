#!/usr/bin/env perl

##########################################
# Check the cluster files
#
# WangMing wangm08(AT)yeah.net
# 2010-07-27 v1.0
##########################################

use strict;
use warnings;
use Getopt::Long;

my ($help,$input);
GetOptions(
	"help"=>\$help,	
	"input=s"=>\$input,
);

if (!defined $input || defined $help) {
	&usage();
	exit 1;
}

open CLU, $input or die "Cannot open file $input.";
open OUT, ">",$input.'.stat.xls' or die "Cannot open file $input.";

my (%h1,$t1,$t2,$t3);
#my @grs = ['cluster1'];

while(<CLU>){
	chomp;
	my @ps = split /\t/;
	$h1{$ps[0]}{'term'} .= '\t'.$ps[0];
	$h1{$ps[0]}{'exp'} += $ps[1];
	my $len = length ($ps[2]);
	$h1{$ps[0]}{'bp'} += $len;	
$t1 ++;
}

foreach my $p (sort keys %h1)	{
    my @h1 = split /\t/, $h1{$p}{'term'};
    my $num = @h1;
	print OUT "$p\t$num\t$h1{$p}{'term'}\t$h1{$p}{'exp'}\t$h1{$p}{'bp'}\n";
#	print OUT "$p\t$h1{$p}\n";
$t2 ++;	
	}
print "$t1\n$t2\n";

close CLU;
close OUT;

sub usage{
	print << "USAGE";
usage: perl cluster_stat.pl -input cluster.txt	
USAGE
}
