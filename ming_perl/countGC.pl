#!/usr/bin/env perl

########################################
# Count GC% of input fasta file
########################################

use strict;
use warnings;

open IN, $ARGV[0] or die "Cannot open file $!.";

my ($len, $GC);

my $count_A = 0;
my $count_T = 0;
my $count_C = 0;
my $count_G = 0;
my $count_N = 0;
my $gap;
<IN>;
while (<IN>){
#     if ($_ =~ /^\>/) {next;}
#       else {
	chomp;
#	$len += length ($_);
        $count_A += ($_ =~ tr/Aa/AA/);
        $count_T += ($_ =~ tr/Tt/TT/);
        $count_C += ($_ =~ tr/Cc/CC/);
        $count_G += ($_ =~ tr/Gg/GG/);
	$count_N += ($_ =~ tr/Nn/NN/);
#           }
}
$len = $count_A + $count_T + $count_C + $count_G + $count_N;
#$GC = ($count_C + $count_G) / ($count_A + $count_T + $count_C + $count_G);
$GC = ($count_C + $count_G) / $len;
$gap = $count_N/$len;

printf "GC%%: %5.2f. \n", $GC*100;

print "A\:$count_A\nT\:$count_T\nC\:$count_C\nG\:$count_G\nN\:$count_N\nGap\:$gap\nLength\:$len\n";

close IN;
