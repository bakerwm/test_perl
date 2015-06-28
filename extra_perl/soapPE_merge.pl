#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Std;

sub help{
    print STDERR <<EOF
Usage: perl  soapPE_merge.pl  -a ref_genome.fa  -i soapPE.txt  -o soapPE_merge.txt  -c  clean.fa  -d soapPE.delete

Note:
    a       : The reference genome. FASTA format.
    i       : The PE soap file.
    o       : The merged PE soap file.
    c       : The merged clean reads in FASTA forma.
Description:
    Merge pair-reads from PE soap file into one read.

EOF
}

my %options = ();
getopts("a:i:o:c:d:", \%options);
if(!defined $options{a} || !defined $options{i} || !defined $options{o} || !defined $options{c} || !defined $options{d}){
    &help();
    exit(1);
}

# Read reference genome.
my $genome;
open F,$options{a} or die;
foreach my $i(<F>){
    next if($i =~ />/);
    chomp($i);
    $genome .= $i;
}
close F;

# Read PE soap file.
open F,$options{i} or die;
open FA,"> $options{c}" or die;
open OUT,"> $options{o}" or die;
open DEL,"> $options{d}" or die;
while(<F>){
    chomp;
    my $t1 = $_;
    chomp(my $t2 = <F>);
    my @a  = &READ($t1);    # id, length, strand, begin.
    my @b  = &READ($t2);
# Check the pair reads were uniq.
    &error($.) if($a[0] ne $b[0]); # The id for read-1 and read-2 should be identical.
    my @m  = &MERGE(@a, @b);    # length, begin, strand.
# Exclude abnomal reads.
    if($m[0] > 300 || $m[0] < 30){
        my $del = join"\t",($m[0],@a, @b);
        print DEL $del,"\n";
        next;
    }
# Print out reads.
    my $seq = substr($genome, ($m[1]-1), $m[0]);
    my $ascii = 'h'x$m[0];
    my @to    = split(/\t/,$t2);
    ($to[0], $to[1], $to[2], $to[5], $to[8]) = ($b[0], $seq, $ascii, $m[0], $m[1]);
    my $out   = join"\t",@to;
    print OUT $out,"\n";
    if($m[2] eq '-'){
        $seq = reverse($seq);
        $seq =~ tr/atcgnATCGN/tagcnTAGCN/;
    }
    print FA ">$a[0]\n$seq\n";
}
close F;
close FA;
close OUT;

sub error{
    my $ln = shift;
    print "Check the line $ln in file \"$options{i}\"\n...It's not a PE soap file!\n";
    exit(1);
}

sub READ{
    my $in = shift;
    my @r  = split(/\t/,$in);
    $r[0]  =~ s/\/\d//g;
    return($r[0], $r[5], $r[6], $r[8]); # id, length, strand, begin
}

sub MERGE{
    my @p = @_; # Read-1 (id, len, str, beg), Read-2 (id, len, str, beg)
    my ($length, $begin, $strand);
    $strand = $p[6];    # Read-2
    my ($read1_beg, $read2_beg)=($p[3],$p[7]);
    my $read1_end = $p[3]+$p[1]-1;
    my $read2_end = $p[7]+$p[5]-1;
    my $NewLength = $read2_end - $read1_beg + 1;
    my ($Newbeg,$Newend) = ($read1_beg,$read2_end);
    if($strand eq '-'){
        $NewLength = $read1_end-$read2_beg+1;
        ($Newbeg,$Newend) = ($read2_beg, $read1_end);
    }
    return($NewLength, $Newbeg,$strand);
}
    
