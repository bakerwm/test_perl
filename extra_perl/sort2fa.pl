#!/user/bin/perl -w
use warnings;
use strict;
use Getopt::Std;

sub help{
    print STDERR <<EOF
Usage: perl  sort2fa.pl  -a  genome.fa  -i input.txt  -o input.fa

Noet:
    a       : genome FASTA file. (for one chromesome).
    i       : Input seq position. (at least 6 lines)
              <ID*> <..> <..> <Begin*> <End*> <Strand*>
    o       : output Sequence file.

EOF
}

my %opt = ();
getopts("a:i:o:", \%opt);
if(!defined $opt{a} || !defined $opt{i} || !defined $opt{o}){
    &help();
    exit(1);
}

open IN,$opt{a} or die;
my $genome;
foreach my $i(<IN>){
    next if($i =~ />/);
    chomp($i);
    $i =~ tr/atcgn/ATCGN/;
    $genome .= $i;
}
close IN;

open IN,$opt{i} or die;
open OUT,"> $opt{o}" or die;
while(<IN>){
    my @tabs = split(/\t/,$_);
    if(@tabs < 6){
        print "File \"$opt{i}\" contain at least 6 lines.\n";
    }

    if($tabs[5] =~ /\w/){
        next if($. == 1);
        print "The 6th part of \"$opt{i}\" should be +/- in line $..\n";
        exit(1);
    }
    my $len = $tabs[4] - $tabs[3] + 1;
    my $begin = $tabs[3] - 1;
    my $fa = substr($genome,$begin,$len);
    my $fa_rc = $fa;
    $fa_rc =~ tr/ATCG/TAGC/;
    $fa_rc = reverse($fa_rc);
    my $seq = ($tabs[5] =~ /\+/)?$fa:$fa_rc;
    my $seq_out = &cut($seq);
    print OUT ">$tabs[0]\n$seq_out";
}
close IN;
close OUT;

sub cut{
    my $seq = shift;
    my $out;
#    for(my $i=0;$i<length($seq);$i+=70){
#        my $line = substr($seq,$i,70);
#        $out .= $line."\n";
#    }
    $out = $seq."\n";
    return $out;
}
