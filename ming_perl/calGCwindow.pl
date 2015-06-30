#!/usr/bin/env perl 
use strict;
use warnings;

sub usage {
    die(qq/
Usage: calGCwindow.pl <window> <in.fa>

Example: calGCwindow.pl 500 Rv.fa
\n/);
}

my $window = shift or usage();
my $infile = shift or usage();

file2window($infile, $window);

sub file2window {
    my ($in, $window) = (@_);
    my $ref = '';
    open my $fh_fa, "< $in" or die "Cannot open file $in, $!\n";
    while(<$fh_fa>) {
        chomp;
        next if(/^>/);
        $ref .= $_;
    }
    close $fh_fa;
    my $ref_length = length($ref);
    my $ref_gc     = calGC($ref);
    for(my $i = 1; $i <= $ref_length; $i += $window) {
        my $start = $i;
        my $end   = $i + $window - 1;
        $end = ($end < $ref_length)?$end:$ref_length;
#        $end < $ref_length || $end = $ref_length;
        my $i_gc  = pos2gc($ref, $start, $end);
        my $out_gc = sprintf"%.4f", ($i_gc - $ref_gc);
        print join("\t", 'chr1', $i, $end, $out_gc) . "\n";
    }
}

######
sub pos2gc {
    my ($seq, $start, $end) = @_;
    my $length = $end - $start + 1;
    $start --; # perl is 0-left most index
    my $sg = substr($seq, $start, $length);
    my $gc = calGC($sg);
    return $gc;
}

sub calGC{
    my $seq = shift(@_);
    die("Find no ATCGN in the seq:\n[$seq]\n") unless ($seq =~ /^[ATCGN]+$/s);
    my $a = $seq =~ tr/Aa/Aa/;
    my $t = $seq =~ tr/Tt/Tt/;
    my $c = $seq =~ tr/Cc/Cc/;
    my $g = $seq =~ tr/Gg/Gg/;
    my $p = sprintf "%.4f", ($g + $c)/($a + $t + $g + $c);
    return $p;
}

