#!/usr/bin/perl  -w
use strict;
use Getopt::Long;
use File::Basename;

my($infile,
   $cutoff,
   $help,
   );

&GetOptions("i=s"  => \$infile,
            "c=s"  => \$cutoff,
            "h"    => \$help,
            );

if( $help or not $infile){
    &help();
    exit(1);
}

sub help {
     print STDERR <<EOF;
: Statistics reads coverage and 5' or 3' ends coverage.
 
Usage: $0  -i match_genome.45SE.txt 
    Options
        -i <file>       : input match_genome file <file>
        -c              : specify cutoff for each read
        -h              : show this help
EOF
}

$infile =~ /genome\.(\w+)\.txt/;
my $id=$1;
my (%hash_5end, 
    %hash_3end,
    %hash_Pos,
    );

$cutoff = ($cutoff)?$cutoff:0;

open (IN,"<$infile") or die;
while(<IN>){
	chomp;
	my @ps=split (/\s+/);
    next unless($ps[6] >= $cutoff);

    my ($end5, $end3) = ($ps[4] eq "+")?($ps[2], $ps[3]):($ps[3], $ps[2]);
	$hash_5end{$ps[4]}->{$end5} += $ps[6];
	$hash_3end{$ps[4]}->{$end3} += $ps[6];
    		
	for(my $k=$ps[2]; $k<=$ps[3];$k++){
		$hash_Pos{$ps[4]}->{$k} += $ps[6];
		}
    		
}
close IN;

my $dir = dirname $infile;
open (END5, ">$dir\/$id\_5end.mapping.txt") or die;
open (END3, ">$dir\/$id\_3end.mapping.txt") or die;
open (PosCov, ">$dir\/$id\_PosCov.mapping.txt") or die;

foreach my $strand (keys %hash_5end){
    foreach my $pos (sort{$a<=>$b} keys %{$hash_5end{$strand}}){
        print END5 "$pos\t$hash_5end{$strand}->{$pos}\t$strand\n";
    }
}
	
foreach my $strand (keys %hash_3end){
	foreach my $pos (sort{$a<=>$b} keys %{$hash_3end{$strand}}){
   		print END3 "$pos\t$hash_3end{$strand}->{$pos}\t$strand\n";
	}
}
	
foreach my $strand(keys %hash_Pos){
	foreach my $pos(sort{$a<=>$b} keys %{$hash_Pos{$strand}}){
   		print PosCov "$pos\t$hash_Pos{$strand}->{$pos}\t$strand\n";
	}
}

close END5;
close END3;
close PosCov;
