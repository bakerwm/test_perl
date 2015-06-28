#!/usr/bin/perl  -w

$usage = "perl  sort2temp.pl  infile  \n";

$f = shift or die $usage;

open F,$f or die;
open OUT,">$f\.temp" or die;
while(<F>){
    chomp;
    next if(/Description/i);
    @a   = split /\t/;
    $length = $a[4] - $a[3] + 1;
    $out = join"\t",($a[0], $a[5], "$a[3]\:Cov", "max\:$a[1]", "$a[4]\:Cov", $length);

    print OUT $out,"\n";
}
close F;
close OUT;
