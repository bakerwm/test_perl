#!/usr/bin/perl -w
use warnings;
use strict;


my $in = shift or die "Need input unique.fa" ;

my $out = $in; $out =~ s/\.fa/.all.fa/;
unlink $out if -e $out;

open F, $in or die;
open OUT, ">>$out" or die;

my $k = 0;
while(<F>){
    next unless(/^\>/);
    my ($id, $count) = /\>(\w+)\s+(\d+)/;
    chomp(my $seq = <F>);
    for(my $i=1;$i<=$count;$i++){
        my $tag = sprintf"%06d",$i;
        my $new_id = $id."\_$tag";
        print OUT ">$new_id\n$seq\n";
        $k ++;
    }
}
close F;
close OUT;

print "output reads: ", $k,"\n";
