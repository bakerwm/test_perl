#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;

my $in = shift or die "Input a fasta file:\n";
check_fasta($in);

sub check_fasta {
    my $seq = shift(@_);
    my $seq_tmp = $seq. ".tmp";
    my $fa_out_num = 0;
    my $seqin = Bio::SeqIO->new(-file    => "< $seq",
                                -format  => 'fasta');
    my $seqout = Bio::SeqIO->new(-file   => "> $seq_tmp",
                                -format  => 'fasta');
    while(my $s = $seqin->next_seq) {
        my $id    = $s->id;
        my $newid = (split /\|/, $id)[3];
        $s->id($newid);
        $s->desc("");
        if(length($s->seq) > 10 ) {
            $seqout->write_seq($s);
            $fa_out_num ++;
        }
    }
    return $fa_out_num;
}
