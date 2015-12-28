#!/usr/bin/perl -w
use warnings;
use strict;

#Description
#
#This script is designed to convert unique fasta file to total.

sub usage {
    die("Usage: fasta_uniq2all.pl <in.fa>\n");
}

usage() if(@ARGV == 0 && -t STDIN);

my $header = "";
my $new_id = "";
my $sequence = "";
while(<>) {
    chomp;
    if($_ =~ /^>(.+)/) {
        if($header ne "") {
             print ( rep_fa($header, $sequence));
       }
        $header   = $1;
        $sequence = "";
    }else {
        $sequence .= $_;
    }
}

print (rep_fa($header, $sequence));

sub rep_fa {
    my $header = $_[0];
    my $seq    = $_[1];
    my $new_seq = "";
    my ($id, $num) = (split /\s+/, $header)[0, 1];
    die("[$header] input fasta shoule be [>ID NUM] format\n") if( ! defined $num || ! $num =~ /%\d+$/);
    for(my $i = 1; $i <= $num; $i ++) {
        my $new_id = sprintf "$id\_%06d", $i;
        $new_seq .= ">" . $new_id . "\n" . $seq . "\n";
    }
    return($new_seq);
}

### END OF FILE ###
