#!/usr/bin/env perl 
### sort (12-col) + count (8-col)
use strict;
use warnings;
use Getopt::Std;

sub usage {
    die("
Usage: pick_TPM.pl [options] <in.txt|STDIN>

Options: -l  <INT>  The number of columns in front of count/TPM start; [12]
         -n  <INT>  The number of libraries; [4]
         -c  <STR>  Calculate the expression of output, 'mean' or 'select'; ['select']
\n");
}

my %opts = (l => 12, n => 4, c => 'select');
getopts("l:n:c:", \%opts);

usage() if(@ARGV == 0 && -t STDIN);

## supposed width of file
my $file_width = $opts{l} + 2 * $opts{n};

while(<>) {
    chomp;
    my @tabs = split /\t/;
    die("Line-$. less than $file_width\n") if(@tabs < $file_width);

    if($opts{c} eq 'mean') {
        mean_exp($_);
    }
    if($opts{c} eq 'select') {
        sel_exp($_);
    }
}

sub sel_exp {
    my $line = $_[0];
### n = 3
    my @tabs = split /\t/, $line;
    my ($count, $tpm);
    if($tabs[2] < 40 ) {
        ($count, $tpm) = @tabs[($opts{l} + 0), ($opts{l} + 1)];
    }elsif($tabs[2] < 80) {
        ($count, $tpm) = @tabs[($opts{l} + 2), ($opts{l} + 3)];
    }elsif($tabs[2] < 140) {
        ($count, $tpm) = @tabs[($opts{l} + 4), ($opts{l} + 5)];
    }else {
        ($count, $tpm) = @tabs[($opts{l} + 6), ($opts{l} + 7)];
    }
    print join("\t", @tabs[0..($opts{l} - 1)], $count, $tpm). "\n";
}

sub mean_exp {
    my $line = $_[0];
### n=3
    my @tabs = split /\t/, $line;
    my $count = sprintf"%.2f", ($tabs[$opts{l} + 0] + $tabs[$opts{l} + 2] + $tabs[$opts{l} + 4]) / $opts{n};
    my $tpm   = sprintf"%.2f", ($tabs[$opts{l} + 1] + $tabs[$opts{l} + 3] + $tabs[$opts{l} + 5]) / $opts{n};
    print join("\t", @tabs[0..($opts{l} - 1)], $count, $tpm) . "\n";
}





### END OF FILE ###
