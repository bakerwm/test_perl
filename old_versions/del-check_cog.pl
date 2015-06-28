### !---------------------------------------
### discard

#!/usr/bin/perl -w
use warnings;
use strict;
use Data::Dumper;

my $usage = "Usage: perl  check_cog.pl  Rv_cog.txt  in_mRNA.txt\n";
my $cog = shift || die "Input COG list:\n$usage\n";
my $in  = shift || die "Input gene list:\n$usage\n";
my $tab = shift;

$tab = 1 if(!$tab);
open F,$in or die;
my %G;
foreach my $i(<F>){
    my @tabs = split(/\t/,$i);
    my $i = $tab - 1;
    $G{$tabs[$i]} = 1;
}
close F;

open F, $cog or die;
foreach my $i(<F>){
    chomp($i);
    my @tabs = split(/\t/,$i);
    my @ps   = split(/\;/,$tabs[2]);
    my @hit = ();
    my $num = 0;
    foreach my $gene (@ps){
        next if(!exists $G{$gene});
        push @hit, $gene;
        $num ++;
    }
    my $gene = join"\,",(@hit);
    my $out  = join"\t",($tabs[0], $tabs[1], $num, $gene);
    print $out,"\n";
}
close F;
