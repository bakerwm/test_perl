### !---------------------------------
### replaced by: chk_seq2rnaz.pl

#!/usr/local/bin/perl
use strict;
use warnings;

use File::Spec::Functions qw(catdir catfile);
use File::Basename qw(dirname basename);
use File::Path qw(make_path remove_tree);
use Getopt::Std;
use POSIX qw(strftime);

use lib qw(/home/wangming/local/lib/perl5); # bioperl module
use Bio::SearchIO;

my %opts = ();
getopts("d:o:", \%opts);
die ("Usage: perl blast2RNAz.pl [-d] <db.fa> [-o] <out.dir> infile.fa\n") unless(@ARGV == 1);
my $db = $opts{d};
my $outdir = $opts{o};
my $infile = shift;

##======Start work dir======##
# 1. Clear output path
die "Warn: remove all files in: [$outdir], and execute again.\n" if(defined (my $check = <$outdir/*>) );
mkdir $outdir unless -d $outdir;
my $fa_path = catdir($outdir, "SeqFA");
make_path($fa_path) unless -d $fa_path;

# 2.Tools
my $blastall = 'blastall'; 
my $formatdb = 'formatdb'; 
my $clustalw2 = 'clustalw2'; 

# 3. Prepare blastdb
die "File not found: $db, $_" unless -e $db;
system "$formatdb -p F -i $db " unless -e "$db\.nhr";
##======END work dir======##

##======Start main script======##
my $t = get_date();
# 1. Run blast (-m7) xml format
print "1. [$t] Run blast\n";
die "File not found: [$infile] \n" unless -e $infile;
my $blast_out = catfile($outdir, "blast.out");
system "$blastall -p blastn -i $infile -d $db -e 1e-2 -m 7 -o $blast_out";

# 2. Parsing blast output file
$t = get_date();
print "2. [$t] Parsing blast output\n";
blast_parser($blast_out);

# 3. Run ClustalW2 and RNAz
$t = get_date();
print "3. [$t] Run clustalW2 and RNAz\n";
my @check_fas = <$fa_path\/*fa>;
die ("Find no fasta files: [$fa_path]") if(@check_fas < 1);
run_clustalw2_rnaz();

# 4. Finish
$t = get_date();
print "4. [$t] Finsh!\n";

##======END main script======##

#####################
#+++ Subroutines +++#
#####################
sub get_date{
    my $datestring = strftime "%Y-%m-%d %H:%M:%S", localtime;
    return $datestring;
}

sub blast_parser {
    my $blast_out = shift(@_); # input the blast output in xml format
    my $in = Bio::SearchIO->new(-format => 'blastxml',
                                -file => $blast_out);
    my $best_hits_file = catfile($outdir, 'best_hits.txt'); ## output best hits txt
    open H, "> $best_hits_file" or die "$!";
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
            $strain = 'Msemg' if ($strain eq 'str');
            my $s_ID = $strain.'_'.$q_name;
            ## output fasta files
            my $seq_fa = catfile($fa_path, "$q_name\.fa");
            open FA, ">> $seq_fa" or die "$!" ;
            print FA "\>$s_ID\n$hsp_string\n";
            close FA;
            my $q_info = join "\t", ($q_name, $s_name, $q_length, $s_length, $frac_identical,
                                     $q_start, $q_end, $s_start, $s_end, $s_evalue, $s_bits);
            print H $q_info, "\n";
        }
    }
    close H;
}

sub run_clustalw2_rnaz{
    my $best_rnaz_bed = catfile($outdir, "best_RNAz.bed");
    open Z, "> $best_rnaz_bed" or die "$!";
    while (my $fa = <$fa_path\/*fa>) {
        sort_fa($fa);
        my $clu_path = catfile($fa_path, basename($fa));
        $clu_path =~ s/\.fa//;
        mkdir $clu_path unless -d $clu_path;
        my $fa_in_path = catfile($clu_path, basename($fa));
        my $clu_log = catfile($outdir, "clustalw2.log");
        rename $fa, $fa_in_path;
        # run clustalw2
        system "$clustalw2 $fa_in_path >> $clu_log";
        # run RNAz
        my $rnaz_path = $clu_path;
        my $clu_aln = my $rnaz_out = $fa_in_path;
        $clu_aln =~ s/\.fa/.aln/;
        $rnaz_out =~ s/\.fa/.bed/;
        my $run_rnaz =  "rnazWindow.pl --min-seqs=2 --max-seqs=6  $clu_aln |\ 
                         RNAz --forward --no-shuffle --cutoff=0.5 | \
                         rnazCluster.pl | rnazIndex.pl --bed | rnazBEDsort.pl > $rnaz_out";
        system "$run_rnaz";
        # filter rnaz_output
        my $best_rnaz = rnaz_parser($rnaz_out);
        next if ($best_rnaz eq '-');
        print Z $best_rnaz, "\n";
    }
    close Z;
}

# move H37Rv to the top
sub sort_fa{
    my $fa = shift(@_); # Input fa file
    my %Seq = ();
    my $fa_out = '';
    my $Rv_id = '';
    open F, $fa or die "$!";
    while (my $in = <F>) {
        chomp($in);
        next unless ($in =~ /^\>/);
        $in =~ s/^\>//g;
        chomp(my $seq = <F> );
        $Seq{$in} = $seq;
        $Rv_id = $in if($in =~ /H37Rv/);
    }
    close F;
    open OUT, "> $fa" or die "$!";
    # Write H37Rv seq
    print OUT "\>$Rv_id\n$Seq{$Rv_id}\n";
    foreach my $n (sort keys %Seq){
        next if ($n eq $Rv_id); # H37Rv
        print OUT "\>$n\n$Seq{$n}\n";
    }
    close OUT
}

# Parsing the BED output of RNAz, get RNAz_score for each input
# * Top RNAz-score fragment
sub rnaz_parser{
    my $in = shift(@_); # input *.bed file
    my %Info = ();
    open F, "<$in" or die "$!";
    while(<F>) {
        chomp;
        next if ($_ eq '');
        my ($id) = (split /\t/)[0]; # select the first hit.
        push @{$Info{$id}}, $_;
    }
    close F;
    my @best_out;
    if((keys %Info) < 1) {
        return '-';  # zero output in bed
    }else {
        my $best_score = 0.5;
        my $best_info = '-';
        foreach my $i (keys %Info){
            foreach my $s (@{$Info{$i}}){
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

blast2RNAz.pl - evaluate the RNAz of input seq using blastn and RNAz.

=head1 SYNOPSIS

perl blastParser.pl -i input.fa -d database.fa  -o result

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
