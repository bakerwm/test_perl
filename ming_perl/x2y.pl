#!/usr/bin/perl -w
use Bio::SeqIO;

$format = 'The following formats supported: fasta, genbank, pir, embl, raw, ace, bsml, swiss';

$usage          = "perl  x2y.pl infile  infileformat  outfile  outfileformat\n\n$format\n";
$infile         = shift or die $usage;
$infileformat   = shift or die $usage;
$outfile        = shift or die $usage;
$outfileformat  = shift or die $usage;

$seq_in=Bio::SeqIO->new(
                        -format => $infileformat,
                        -file   => "<$infile",
                        );

$seq_out=Bio::SeqIO->new(
                        -format => $outfileformat,
                        -file   => ">$outfile",
                        );

while( $inseq = $seq_in->next_seq){
    $seq_out->write_seq($inseq);
}
