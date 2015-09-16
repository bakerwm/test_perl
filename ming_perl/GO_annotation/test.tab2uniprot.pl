#!/usr/bin/env perl

use strict;
use warnings;

die("Usage: tab2uniprot.pl in.tab\n") if(@ARGV == 0 && -t STDIN);

while(<>) {
    chomp;
    my @tabs = split /\t/;
    my @genes = split /\s+/, $tabs[4];
    for my $g (@genes) {
        print join("\t", $g, @tabs) . "\n";
    }
}
