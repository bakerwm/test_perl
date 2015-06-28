### !-------------------------------------------------
### delete, replaced by: bedtools

#!/usr/bin/perl  -w
use warnings;
use strict;
use Getopt::Std;
use File::Basename;

sub help{
    print STDERR <<EOF
Usage: perl  ExtFlank.pl  -g genome.fna  -i input.txt  -t 5end -o out.dir

Note:
Extract the 5/3 flanking sequences of input txt.
    g       : The genome sequence file. [FASTA]
    i       : The sorted file; Tab-delimited.
    t       : Type, [5end/3end/both/].
              default [5end]
    o       : Directory of output file. 
              default [./]
    s       : The length of sequence within input sequence. [integer]
              default [0], [A] indicate the whole input sequence.
    e       : The length of sequence upstream/downstream of input sequence. [integer]
              default [100].
EOF
}

my %options = ();
getopt("g:i:t:o:s:e:", \%options);
if(!defined $options{g} || !defined $options{i}){
    &help();
    exit(1);
}
$options{t} = (defined $options{t})?$options{t}:'5end';
$options{o} = (defined $options{o})?$options{o}:'./';
$options{s} = (defined $options{s})?$options{s}:0;
$options{e} = (defined $options{e})?$options{e}:100;

unless($options{t} eq '5end' || $options{t} eq '3end' || $options{t} eq 'both'){
    print "Check parameter: -t, which should be [5end/3end/both].\nDefault is: 5end\n";
    exit(1);
}
 
#######################
##  Read genome seq  ##
#######################
open F,$options{g} or die "cannot open file: $options{g}, $! \n";
my $RefSeq;
my $RefLength;
while(<F>){
    chomp;
    next if(/\>/);
    $RefSeq .= $_;
    $RefLength += length($_);
}
close F;

open F, $options{i}  or die "cannot open file: $options{i}, $! \n";
my $printPromSeq;
my $printTermSeq;
my $printPromTxt;
my $printTermTxt;
while(<F>){
    chomp;
    my @ps = split(/\t/,$_);
    my $SeqLength=$ps[4]-$ps[3]+1;
    my $InLength=($options{s} eq 'A')?$SeqLength:$options{s};
    my $OutLength=$options{e};
    my $LeftStart   = $ps[3] - $options{e};
       $LeftStart   = 0 if($LeftStart < 0);
    my $LeftEnd     = $ps[3] + $InLength - 1;
    my $RightStart  = $ps[4] - $InLength + 1;
    my $RightEnd    = $ps[4] + $options{e};
    my $CutLength   = $InLength + $OutLength;
    my $LeftSeq     = &Seq($LeftStart, $CutLength, $ps[5]);
    my $RightSeq    = &Seq($RightStart, $CutLength, $ps[5]);
    my ($PromSeq, $TermSeq, $PromTxt, $TermTxt);
    ($PromSeq, $TermSeq) = ($LeftSeq, $RightSeq);
    $PromTxt = join"\t",("$ps[0]\_prom", 'prom', $CutLength, $LeftStart, $LeftEnd, $ps[5], "\n");
    $TermTxt = join"\t",("$ps[0]\_term", 'term', $CutLength, $RightStart, $RightEnd, $ps[5], "\n");
    if($ps[5] eq '-'){
        ($PromSeq, $TermSeq) = ($RightSeq, $LeftSeq);
        $PromTxt = join"\t",("$ps[0]\_prom", 'prom', $CutLength, $RightStart, $RightEnd, $ps[5], "\n");
        $TermTxt = join"\t",("$ps[0]\_term", 'term', $CutLength, $LeftStart, $LeftEnd, $ps[5], "\n");
    }
    my $PromLine   = '>'.$ps[0].'_prom'."\n".$PromSeq."\n";
    my $TermLine   = '>'.$ps[0].'_term'."\n".$TermSeq."\n";
# 
    ($PromLine, $PromTxt) = ('', '') if(length($PromSeq) < $options{e} || length($PromSeq) > $CutLength);
    ($TermLine, $TermTxt) = ('', '') if(length($TermSeq) < $options{e} || length($TermSeq) > $CutLength);
    $printPromSeq .= $PromLine;
    $printTermSeq .= $TermLine;
    $printPromTxt .= $PromTxt;
    $printTermTxt .= $TermTxt;
}
close F;

my $filename = basename($options{i});
   $filename =~ s/\.txt//g;
my ($flank5, $flank3) = ("$filename\_5flank.fa", "$filename\_3flank.fa");
my ($Txt5, $Txt3) = ("$filename\_5flank.txt", "$filename\_3flank.txt");

if($options{t} eq 'both'){
    open PRO, "> $flank5" or die "cannot open file: $flank5, $! \n";
    open TER, "> $flank3" or die "cannot open file: $flank3, $! \n";
    open T1,"> $Txt5" or die;
    open T2,"> $Txt3" or die;
    print PRO $printPromSeq;
    print TER $printTermSeq;
    print T1  $printPromTxt;
    print T2  $printTermTxt;
    close PRO;
    close TER;
    close T1;
    close T2;
}elsif($options{t} eq '3end'){
    open TER, "> $flank3" or die "cannot open file: $flank3, $! \n";
    open T2,"> $Txt3" or die;
    print TER $printTermSeq;
    print T2  $printTermTxt;
    close TER;
    close T2;
}else{
    open PRO, "> $flank5" or die "cannot open file: $flank5, $! \n";
    open T1,"> $Txt5" or die;
    print PRO $printPromSeq;
    print T1  $printPromTxt;
    close PRO;
    close T1;
}

sub Seq{
    my ($start, $length, $strand) = @_;
    my $seq = substr($RefSeq, $start, $length);
    if($strand eq '-'){
        $seq = reverse($seq);
        $seq =~ tr/ATCGatcg/TAGCtagc/;
    }
    return $seq;
}
