#!/usr/bin/perl  -w
use warnings;
use strict;
use Getopt::Std;

sub help{
    print STDERR <<EOF
Usage: perl  overlap.pl  -i  gene.list  -d database.txt  -r 45  -c 0.5

Note:
    i       : The input list should contain at least 6 lines.
              <ID*> <..> <..> <Begin*> <End*> <Strand*>...
    d       : The database file. at least 6 lines.
              <ID*> <..> <..> <Begin*> <End*> <Strand*>...
    r       : The estimate max length of input sequences.
              Default [1000].
    c       : Cutoff for the overlap to each sequence.
              Default [0.5]
EOF
}

# Set parameters.
my %options = ();
getopt("i:d:r:c:",\%options);
if(!defined $options{i} || !defined $options{d}){
    &help();
    exit(1);
}
$options{r} = ($options{r})?$options{r}:1000;
$options{c} = ($options{c})?$options{c}:0.5;

# Read dataset file.
my %Match;
open IN, $options{d} or die;
foreach my $i(<IN>){
    chomp($i);
    my @tabs = split(/\t/,$i);
    my @ps;
    for(my $k=0;$k<12;$k++){        # keep no more than 12 lines in file.
        last if(!exists $tabs[$k]);
        push @ps,($tabs[$k]);
    }
    my $out = join"\t",@ps;
    $Match{$tabs[3]}->{$out} = 1;   # Begin position as key1.
}
close IN;

# Read input file and output file.
open F, $options{i} or die;
open OUT_S, "> overlap_sense.list" or die;
#open OUT_A, "> overlap_anti.list" or die;
open OS, "> overlap_sense.txt" or die;
#open OA, "> overlap_anti.txt" or die;
foreach my $k (<F>){
    my $hit_sense = my $hit_anti = ">$k";   # as list.
    chomp($k);
    my @ps = split(/\t/, $k);
# check Head_lines;
    if(!$ps[3]=~/^\d+$/ || !$ps[4]=~/^\d+$/ || !$ps[5]=~/\+|\-/){ # begin, end, strand.
        print "line",$.,"\:",$k,"\n"; # 
        next;
    }
    my $tag = (defined $ps[12])?$ps[12]:'sRNA';
    my ($ge_begin, $ge_end, $ge_strand) = ($ps[3], $ps[4], $ps[5]);
    my $ge_len = $ge_end - $ge_begin + 1;
    my $ge_pos = join"\:",($ge_len, $ge_begin, $ge_end, $ge_strand);
# Search the dataset.
    for(my $j=($ps[3]-$options{r}); $j<=$ps[4]; $j++){
        next unless(exists $Match{$j});
        my @reads = keys(%{$Match{$j}});
        foreach my $m(@reads){
            my @r = split(/\t/,$m);
# The same strand.
            next unless($r[5] eq $ge_strand);
# Have overlap.
            next unless($r[3]<=$ge_end && $r[4]>=$ge_begin);

            my $r_len = $r[4] - $r[3] + 1;
            my $r_pos = join"\:",($r_len, $r[3], $r[4], $r[5]);
            my $overlap_a = $r[4] - $ge_begin + 1;
            my $overlap_b = $ge_end - $r[3] + 1;
            my $gap = ($r[3]<$ge_begin)?(($overlap_a>$ge_len)?$ge_len:$overlap_a):(($overlap_b>$r_len)?$r_len:$overlap_b);
            my $ge_percent = sprintf "%.2f",($gap/$ge_len);
            my $r_percent  = sprintf "%.2f",($gap/$r_len);
# Check: more than 1/2 of either sequence.
            next if($ge_percent < $options{c} && $r_percent < $options{c});
# Check: Stat sense and antisense strand.
            my $overlap_out = join"\t",($ps[0], $ge_pos, $r[0], $r_pos, $gap, $ge_percent, $r_percent);
#            if($ge_strand eq $r[5]){
            print OS $overlap_out,"\n";
            $hit_sense .= "$m\t$ps[0]\:$tag\n";
#           }else{
#                print OA $overlap_out,"\n";
#                $hit_anti .= "$m\t$ps[0]\:$tag\n";
#            }
#            }else{
#                next;
#            }
        }
    }
    print OUT_S $hit_sense;
#    print OUT_A $hit_anti;
}
close F;
close OS;
#close OA;
close OUT_S;
#close OUT_A;
