#! /user/bin/env perl

#######################################
# Convert BAM to wig file
#
#######################################

use strict;
use warnings;

unless(defined @ARGV){
	print "Used for one chromosome results.
usage: perl bam2wig.pl  file.bam  name  description > file.wig\n";
	exit;
	}
#track type=wiggle_0 name=worm description=dfdsfsdvariableStep chrom=chr

my $f=shift;
my $name=shift;
my $des=shift;
open IN, "$f" || die "Cannot open file $f: $!\n";

my $chr;
my $g=<IN>;
chomp $g;
my @ps = split /\t/, $g;
$chr=$ps[0];
print "track type=wiggle_0 name=$name description=$des\n";
print "variableStep chrom=chr\n";
print "$ps[1]\t$ps[2]\n";

while(<IN>){
	chomp;
	my @ps=split /\t/;
	if($ps[0] eq  $chr){
		print "$ps[1]\t$ps[2]\n";
		}
		else{
			  print "variableStep chrom=$ps[0]\n";
			  print "$ps[1]\t$ps[2]\n";
			  $chr=$ps[0];
                   }	
	}
	
close IN;	
