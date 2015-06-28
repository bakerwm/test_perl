#!/usr/bin/perl -w

use strict;
use warnings;

my $pl = 'sort2bed.pl';
my $in_dir = 'example';
my @files = <$in_dir/*>;

die "Files not found at: [$in_dir]" if (@files < 5);

my $genome = 'example/test.fa';
my $feature = 'gene';
my @types = ('gff', 'ptt', 'bed', 'sort', 'fa');

foreach my $f (@files){
    next unless($f =~ /.(gff|ptt|bed|txt)$/);
    my ($in_type) = $f =~ /\.(\w+)$/;
    $in_type = 'sort' if($in_type eq 'txt');
    foreach my $out_type (@types){
        next if($in_type eq $out_type);
        my $t = $in_type.'2'.$out_type;
        my $outFile = $t.'.out';
        print "perl sort2bed.pl -g $genome -f $feature -i $f -t $t -o $outFile \n";
    }
}
