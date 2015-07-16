#!/usr/bin/env perl 

### sort (12-col) + count (8-col)

use strict;
use warnings;


sub usage {
    die("Usage: pick_TPM.pl <in.txt|STDIN>\n");
}

usage() if(@ARGV == 0 && -t STDIN);

while(<>) {
    chomp;
    my @tabs = split /\t/;
    die("Line-$. less than 20-col\n") if(@tabs < 20);
    my ($count, $tpm);
    if($tabs[2] < 40 ) {
        ($count, $tpm) = @tabs[12,13];
    }elsif($tabs[2] < 80) {
        ($count, $tpm) = @tabs[14,15];
    }elsif($tabs[2] < 140) {
        ($count, $tpm) = @tabs[16,17];
    }else {
        ($count, $tpm) = @tabs[18,19];
    }
    print join("\t", @tabs[0..11], $count, $tpm). "\n";
}
