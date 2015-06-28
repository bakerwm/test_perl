### !---------------------------------------------------
### replaced by: rbh_run_blast.pl

#!/usr/bin/env perl

#########################################################
# search_seq.pl
# 
# Designed to search the best hit sequences in database
# using BLAST
#
#########################################################

use strict;
use warnings;
use File::Basename qw(dirname basename);
use File::Spec::Functions qw(catdir catfile);
use Bio::SearchIO;

use Getopt::Std;

my $result = run_blast();

# Subroutines #
sub run_blast{
    my ($in, $db) = get_para();
    # check blastn and makeblastdb in your $PATH
    my $chk_blastn = `which blastn`;
    my $chk_makeblastdb = `which makeblastdb`;
    if($chk_blastn =~ /which\:\sno\sblastn\sin/ || $chk_makeblastdb =~ /which\:\sno\smakeblastdb\sin/){
        die "[blastn|makeblastdb] not found in your PATH. Install blast+ and 
            add it to your PATH. 
            see: [ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/] \n";    
    }
    # make db
    my $node1 = system "formatdb -p F -i $db ";
#    my $node1 = system "makeblastdb -in $db -dbtype nucl -logfile $db.log";
    # run blast
    my $out = "tmp.blastout." . int(rand(10000));
    my $node2 = system "blastall -p blastn -i $in -d $db -m 7 -o $out";
#    my $node2 = system "blastn -query $in -db $db -outfmt 5 -out $out"; # xml output
    my $results = blast_best_hits($out);
    print $results,"\n";
    unlink($out);
}

sub get_para{
    my %opts = ();
    getopts("d:h:", \%opts);
    usage() if(@ARGV == 0);
#    die ("Usage: search_seq.pl [-d] <db.fa> <infile.fa>\n") if (@ARGV == 0);
    my $db = $opts{d};
    die "[$db] not found" unless -e $db;
    my $infile = shift(@ARGV);
    return ($infile, $db);
}

sub blast_best_hits{
    my $infile = shift(@_); # the balstout file
    my $in = Bio::SearchIO->new(-format => 'blastxml',
                                  -file   => $infile);
    my @output = ();
    while(my $result = $in->next_result) {
        my $q_name = $result->query_description; # why not query_name?
        my $q_length = $result->query_length;
        my $count = 0;
        while(my $hit = $result->next_hit) {
            my $hsp_num = $hit->num_hsps; # ?? number of hits
            my $hsp = $hit->next_hsp; # choose the first hsp
            my $s_length = $hit->length;
            my $s_name = $hit->name;
            my $mapped = $hsp->length('hit');
            my $s_evalue = $hsp->evalue;
            my $s_bits = $hsp->bits;
            my $frac_identical = sprintf"%.1f", $hsp->frac_identical * 100;
            my $hsp_string = $hsp->hit_string;

            my ($q_start, $q_end) = ($hsp->start('query'), $hsp->end('query'));
            my ($s_start, $s_end) = ($hsp->start('hit'), $hsp->end('hit'));
            my ($q_strand, $s_strand) = ($hsp->strand('query'), $hsp->strand('hit'));
            ($q_start, $q_end) = ($q_end, $q_start) if($q_strand < 0);
            ($s_start, $s_end) = ($s_end, $s_start) if($s_strand < 0);
            $q_strand = ($q_strand < 0)?'-':'+';
            $s_strand = ($s_strand < 0)?'-':'+';
            my $map_length = $q_end - $q_start + 1;
#            my $line = join "\t", ($q_name, $s_name, $mapped, $s_evalue, $s_bits, $frac_identical);

             my $line = join "\t", ($q_name, $s_name, $frac_identical, $map_length, '-', '-',
                                    $q_start, $q_end, $s_start, $s_end, $s_evalue, $s_bits,
                                    $q_length, $s_length, $q_strand, $s_strand);
            push @output, $line;
            $count ++;
        }
        next if($count == 0);
    }
    unlink($in);
    return join "\n", @output;
}

sub usage {
    die(qq/
Usage: search_seq.pl [options] input.fa

Options: -d     : The database file [FASTA]

input.fa FASTA file
\n/);
}

__END__

=head1 NAME

SearchSeq.pl - Search the best hit in database using blast.

=head1 SYNOPSIS

perl SearchSeq.pl -i in.fa -d db.fa -o out.txt

=head1 OPTIONS
:wq

=over 8

=item B<-i> Input the seq in fasta format

Choose the input file in FASTA format

=item B<-d> databse

Choose the databse file in FASTA format

=item B<-o> outfile

Choose the file to store the output

=item B<-h> help

Print this help

=back

=cut
