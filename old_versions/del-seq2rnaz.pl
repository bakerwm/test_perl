### !---------------------------------------
### replaced by: chk_seq2rnaz.pl

#!/usr/bin/perl -w
use strict;
use warnings;

#############################
# RNAz evaluation for seqs
# Using Mtb-6 genomes as ref
# 2015-05-05
#############################

use File::Spec::Functions qw(catdir catfile);
use File::Basename qw(dirname basename);
use File::Which;
use File::Path qw(make_path remove_tree);
use Getopt::Std;
use POSIX qw(strftime);

use lib qw(/home/wangming/local/lib/perl5); # bioperl module
use Bio::SearchIO;


my $db;
my $input_fa;
my $outpath;
my $blastall = 'blastall';
my $formatdb = 'formatdb';
my $clustalw2 = 'clustalw2';

parse_para();
my @hit_seqs = exe_blast($input_fa);
my $best_rnaz = exe_clustalw_rnaz(@hit_seqs);

sub parse_para{
    my %opts = ();
    getopt("d:o:", \%opts);
    die ("Usage: perl $0 [-d] <db.fa> [-o] <out.dir> infile.fa\n") if(@ARGV != 1);
    die "[-d] need database fa" unless (-e $opts{d});
    die "[-o] need specify the output dir" unless defined $opts{o};# (-e $opts{o});
    $input_fa = shift(@ARGV);
    die "[$input_fa] file not exists" unless -e $input_fa;
    $db = $opts{d};
    $outpath = $opts{o};
    prep_tools();
    prep_wkpath();
}

sub prep_tools{
#    my $env_path = qx(echo \$PATH);
    foreach my $t ('blastall', 'formatdb', 'clustalw2', 'RNAz'){
        die "[$t] not found in your \$PATH" unless (which($t));
    }
}

sub prep_wkpath{
    my $fa_path = catdir($outpath, "SeqFA");
    make_path($fa_path) unless -d $fa_path;
    # clear previous seq in SeqFA
    while(my $f = <$fa_path/*>){
        remove_tree($f);                         
    }
}

sub exe_blast{
    my $fa = shift(@_);
    # formatdb
    system "$formatdb -p F -i $db";
    my $fa_blast_out = catfile($outpath, "blast.out");
    system "$blastall -p blastn -i $fa -d $db -e 1e-2 -m 7 -o $fa_blast_out";
    my @fa_ids = parse_blast($fa_blast_out);
    return @fa_ids;
}

sub parse_blast{ # extract one best hit for each genome
    my $blast_out = shift(@_); # xml format
    my $in = Bio::SearchIO->new(-format => 'blastxml',
                                -file => $blast_out);
    my $best_hits_file = catfile($outpath, 'best_hits.txt'); ## output best hits txta
    my %seqs;
    open my $fh_hits, "> $best_hits_file" or die "$!";
    while (my $result = $in->next_result) {
        next if($result->num_hits == 0);
        my $q_name = $result->query_name;
        my $q_length = $result->query_length;
        my $q_acc = $result->query_accession;
        my $count = 0;
        while(my $hit = $result->next_hit) {
            my $s_name = $hit->name;
            my $s_length = $hit->length;
            my $s_sig = $hit->significance;
            my $s_description = $hit->description;
            my $hsp_num = $hit->num_hsps;
            # Choose the first HSP hit;
            my $hsp = $hit->next_hsp;
            my ($q_start, $q_end) = ($hsp->start('query'), $hsp->end('query'));
            my ($s_start, $s_end) = ($hsp->start('hit'), $hsp->end('hit'));
            my $s_evalue = $hsp->evalue;
            my $s_bits = $hsp->bits;
            my $frac_identical = sprintf "%.1f", $hsp->frac_identical * 100;
            my $frac_conserved = sprintf "%.1f", $hsp->frac_conserved * 100;
            my $hsp_string = $hsp->hit_string;
            # Create new strain name;
            my ($strain) = $s_description =~ /\w+\s\w+\s(\w+)/;
#            $strain = 'Msmeg' if ($strain eq 'str');
            my $s_ID = $strain.'_'.$q_name;
            ## output fasta files
            $seqs{$q_name} ++;
            my $seq_fa = catfile($outpath, "$q_name\.fa");
            open my $fh_fa, ">> $seq_fa" or die "$!" ;
            print $fh_fa "\>$s_ID\n$hsp_string\n";
            close $fh_fa;
            my $q_info = join "\t", ($q_name, $s_name, $q_length, $s_length, $frac_identical,
                                     $q_start, $q_end, $s_start, $s_end, $s_evalue, $s_bits);
            print $fh_hits $q_info, "\n";
        }
    }
    close $fh_hits;
#    return $best_hits_file;
    return (keys %seqs);
}

sub exe_clustalw_rnaz{
    my @fa_ids = @_;
    my $best_rnaz_bed = catfile($outpath, "best_RNAz.bed");
    open my $fh_rnaz, "> $best_rnaz_bed" or die "$!";
    foreach my $f (@fa_ids) {
        my $fa = catfile($outpath, $f.".fa");
        sort_fa($fa);
        my $fa_path = catdir(dirname($fa), "SeqFA");
        my $clw_path = catdir($fa_path, basename($fa));
        $clw_path =~ s/\.fa//;
        my $fa_new   = catfile($clw_path, basename($fa));
        my $clw_log  = catfile($clw_path, 'clustalw2.log');
        make_path($clw_path) unless -d $clw_path;
        my $rnaz_path = $clw_path;
        rename $fa, $fa_new;
        my $clw_aln  = my $rnaz_out = $fa_new;
        $clw_aln  =~ s/\.fa/\.aln/;
        $rnaz_out =~ s/\.fa/\.bed/;
        # perform ClustalW2
        system "$clustalw2 $fa_new >> $clw_log";
        my $run_rnaz = join (" ", ('rnazWindow.pl --min-seqs=2 --max-seqs=6', $clw_aln, 
                    '| RNAz --forward --no-shuffle --cutoff=0.5', 
                    '| rnazCluster.pl | rnazIndex.pl --bed | rnazBEDsort.pl >', $rnaz_out));
        system "$run_rnaz";
        # filter rnaz_output
        my $best_rnaz = parse_rnaz($rnaz_out);
        next if ($best_rnaz eq '-');
        print $fh_rnaz $best_rnaz, "\n";
    }
    close $fh_rnaz;
    return $best_rnaz_bed;
}

sub show_date{
    my $date = strftime "%Y-%m-%d %H:%M:%S", localtime;
    return $date;
}

# move H37Rv to the first seq
sub sort_fa{
    my $fa = shift(@_); # Input fa file
    my %seq_new = ();
    my $fa_out = '';
    my $Rv_id = '';
    open my $fh_fa, "< $fa" or die "$!";
    while (my $in = <$fh_fa>) {
        chomp($in);
        next unless ($in =~ /^\>/);
        $in =~ s/^\>//g;
        chomp(my $seq = <$fh_fa> );
        $seq_new{$in} = $seq;
        $Rv_id = $in if($in =~ /MG1655/);
#        $Rv_id = $in if($in =~ /H37Rv/);
    }
    close $fh_fa;
    open my $fh_new, "> $fa" or die "$!";
    # Write H37Rv seq
    print $fh_new "\>$Rv_id\n$seq_new{$Rv_id}\n";
    foreach my $n (sort keys %seq_new){
        next if ($n eq $Rv_id); # H37Rv
        print $fh_new "\>$n\n$seq_new{$n}\n";
    }
    close $fh_new;
}

# Parsing the BED output of RNAz, get RNAz_score for each input
# * Top RNAz-score fragment
sub parse_rnaz{
    my $in = shift(@_); # input *.bed file
    my %info = ();
    open my $fh_rnaz, "< $in" or die "$!";
    while(<$fh_rnaz>) {
        chomp;
        next if ($_ eq '');
        my ($id) = (split /\t/)[0]; # select the first hit.
        push @{$info{$id}}, $_;
    }
    close $fh_rnaz;
    my @best_out;
    if((keys %info) < 1) {
        return '-';  # zero output in bed
    }else {
        my $best_score = 0.5;
        my $best_info = '-';
        foreach my $i (keys %info){
            foreach my $s (@{$info{$i}}){
                my ($score) = (split /\t/, $s)[-1];
                ($best_score, $best_info) = ($score, $s) if($score >= $best_score);
            }
            push @best_out, $best_info;
        }
        return join "\n", @best_out;
    }
}

__END__

=head1 NAME

seq2rnaz.pl - evaluate the RNAz of input seq using blastn and RNAz.

=head1 SYNOPSIS

perl blastParser.pl -d database.fa  -o result input.fa 

=head1 OPTIONS

=over 8

=item B<-i>, B<--input>

Choose the input should be in FASTA format.

=item B<-d>, B<--database>

Select the database in FASTA format (NCBI, *.fna format)

=item B<-o>, B<--outdir>

The output directory. defaule [out]

=item B<-h>, B<--help>

Print this help

=back

=head1 DESCRIPTION

B<blast2RNAz.pl>, the analysis consist of three main steps:

1. Run blast: extract the homologe sequences from database using blastn.

2. Run ClustalW2: Using the default para.

3. Run RNAz: Default para.

The output dir contain:

Outdir/   

|-best_hits.txt : The blastn output with best hits.    
    
|-best_RNAz.bed : The RNAz analysis for each input seqs.

|-SeqFA/

|  |- Seq01/ Seq01.fa Seq.aln, Seq.bed, Seq.dnd

|  |- Seq02/ Seq02.fa Seq.aln, ...

|  |- ...

=head1 AUTHOR

Wang Ming, wangmcas@gmail.com

=cut
