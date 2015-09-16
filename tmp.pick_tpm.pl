#!/usr/bin/env perl 

### sort (12-col) + count (8-col)

use strict;
use warnings;
use Getopt::Std;

sub usage {
    die("
Usage: pick_TPM.pl [options] <in.txt|STDIN>

Options: -l  <INT>  The number of columns in front of count/TPM start; [12]
\n");
}

my %opts = (l => 122);
getopts("l:", \%opts);

usage() if(@ARGV == 0 && -t STDIN);

## supposed width of file
my $file_width = $opts{l} + 8;

while(<>) {
    chomp;
    my @tabs = split /\t/;
    die("Line-$. less than $file_width\n") if(@tabs < $file_width);
    my ($count, $tpm);
    if($tabs[2] < 40 ) {
        ($count, $tpm) = @tabs[$opts{l}, ($opts{l} + 1)];
    }elsif($tabs[2] < 80) {
        ($count, $tpm) = @tabs[($opts{l} + 2), ($opts{l} + 3)];
    }elsif($tabs[2] < 140) {
        ($count, $tpm) = @tabs[($opts{l} + 4), ($opts{l} + 5)];
    }else {
        ($count, $tpm) = @tabs[($opts{l} + 6), ($opts{l} + 7)];
    }
    print join("\t", @tabs[0..($opts{l} - 1)], $count, $tpm). "\n";
}
