#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename qw(basename);


my $dir = shift or die "Input the dir of summary files";

my @files = glob"$dir/*summary";

for my $s (@files) {
    my $s_name = basename($s);
    $s_name =~ s/_count.summary//;
    my @counts = parse_summary($s);
    print join("\n", $s_name, @counts) . "\n";
}


sub parse_summary {
    my $in = shift(@_);
    my $line = '';
    open my $fh_in, "< $in" or die "Cannot opne $in, $!\n";
    while(<$fh_in>) {
        $line .= $_;
    }
    close $fh_in;
    my @out = ();
    my @files = split /Status/, $line;
    shift(@files); # trim blank line
    for my $f (@files) {
        my @rs = split /\n/, $f;
        my $header = shift(@rs);
        my $hname  = basename($header);
        my @nums = ();
        for my $r (@rs) {
            my $n = (split /\s+/, $r)[1];
            push @nums, $n;
        }
        push @out, join("\t", $hname, @nums);
    }
    return @out;
}
