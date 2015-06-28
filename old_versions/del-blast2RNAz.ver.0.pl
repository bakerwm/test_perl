### !-------------------------------
### replaced by: chk_seq2rnaz.pl

#!/usr/local/bin/perl
use strict;
use warnings;

use lib qw(/home/wangming/localperl/share/perl5 ); # local module lib

use File::Spec::Functions qw(catdir catfile);
use File::Basename qw(dirname basename);
use File::Which qw(which);
#use Cwd qw(abs_path cwd);
use POSIX qw(strftime);

use lib qw(/home/wangming/local/lib/perl5);

use Bio::SearchIO;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my $infile = '';
my $db = '';
my $outdir = '';
my $help = '';
my $man = '';

GetOptions('i|input=s' => \$infile,
           'd|database=s' => \$db,
           'o|outdir=s' => \$outdir,
           'h|help' => \$help,
           'man' => \$man
           ) or pod2usage(-verbose => 0);

pod2usage(-verbose => 1) if ($help);
pod2usage(-verbose => 2) if ($man);
pod2usage("$0 : [-i|--input $infile] Not found") unless (-e $infile);
pod2usage("$0 : [-d|--database $db] Not found") unless (-e $db);

$outdir = 'out' if ($outdir eq '');
mkdir $outdir unless -d $outdir;

### Prepare result dir
my @tmp_fa = <$outdir\/*.fa>;
if(@tmp_fa > 0){
    foreach my $f (@tmp_fa){
        print "Delete file: $f\n";
        system "rm -f $f";
    }
}

### blast tools
my $blastall = which('blastall');
my $formatdb = which('formatdb');
pod2usage(-verbose => 1, -message => "Cannot find blastall or formatdb in ENV") unless (defined $blastall && defined $formatdb);

### Prepare database
my $db_num = '';
$db_num = `grep '>' $db | wc -l`;  chomp($db_num);

die "File not found: $db, $_" unless -e $db;
system "$formatdb -p F -i $db " unless -e "$db\.nhr";

### Run blast (-m7)
die "File not found: $infile, $_" unless -e $infile;
system "$blastall -p blastn -i $infile -d $db -e 1e-2 -m 7 -o tmp.blst";

### Parsing blast output file
my $in = Bio::SearchIO->new(-format => 'blastxml',
                            -file => 'tmp.blst');
# Write best hits to hit file
# Write fasta file to separate files
my $tab_out = catfile($outdir, 'best_hits.txt');
open TAB, "> $tab_out" or die "$!";
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
### output fasta files
        my $file_out = catfile($outdir, "$q_name\.fa");
        open FA, ">>$file_out" or die ;
        print FA "\>$s_ID\n$hsp_string\n";
        close FA;

### output tab-delimted file
        my $tab_out = join "\t", ($q_name, $s_name, $q_length, $s_length, $frac_identical,
                                  $q_start, $q_end, $s_start, $s_end, $s_evalue, $s_bits);
        print TAB $tab_out, "\n";            
    }
}

### Run clustalW2 and RNAz test for each seq
unlink 'clustalw2.log' if -e 'clustalw2.log';
my $RNAz_out = catfile($outdir, 'best_RNAz.bed');
open RNAz_OUT, "> $RNAz_out" or die;

my $clustalw2 = which('clustalw2');
die "Cannot find clustalw2" unless(defined $clustalw2);
my @fa_files = <$outdir\/*.fa>;
foreach my $fa (@fa_files) {
# Show what is running
    my $datestring = strftime "%Y-%m-%d %H:%M:%S", localtime;
    print "$datestring\tProcessing $fa\n";
    &reorderFA($fa);
    my $file_name = basename($fa);
    $file_name =~ s/\.fa$//;
    my $fa_subdir = catdir($outdir, $file_name);
    mkdir $fa_subdir unless -d $fa_subdir;
    system "mv -f $fa $fa_subdir";
# Run ClustalW2
    system "$clustalw2 $fa_subdir\/*.fa >> clustalw2.log";
#    my $aln_out = catfile($fa_subdir, "$file_name\.aln");
#    system "clustalo -i $fa_subdir\/*.fa --outfmt=clu -o $aln_out --force";
# Run RNAz
    my $in_prefix = catfile($fa_subdir, "$file_name");
    my $in_aln = $in_prefix.'.aln';
    my $aln2window = $in_prefix.'.windows';
    my $RNAzout = $in_prefix.'.rnazout';
    my $RNAzdata = $in_prefix.'.dat';
    my $RNAzbed = $in_prefix.'.bed';
    my $RunRNAz = "rnazWindow.pl --min-seqs=2 --max-seqs=6  $in_aln |\
                   RNAz --forward --no-shuffle --cutoff=0.5 | \
                   rnazCluster.pl | rnazIndex.pl --bed | rnazBEDsort.pl > $RNAzbed";
    system "$RunRNAz";

# Filter the RNAz output (Find the best score)
    my $best_output = '';
    $best_output = &RNAz_Parse($RNAzbed);
    print RNAz_OUT $best_output, "\n";   
}

close RNAz_OUT;

# move H37Rv to the top
sub reorderFA{
    my $fa = shift(@_);
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
sub RNAz_Parse{
    my $bed = shift(@_);
    my %Info = ();
    open F, $bed or die "$!" ;
    while(<F>) {
        chomp;
        next if ($_ eq '');
        my ($id) = (split /\t/)[0];
        push @{$Info{$id}}, $_;
    }
    close F;
    
    my @best_out;
    if((keys %Info) < 1) {
        return '';  # zero output in bed
    }else {
        my $best_score = 0.5;
        my $best_info = '';
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

Outdir    
|    
|--best_hits.txt : The blastn output with best hits.    
|    
|--best_RNAz.bed : The RNAz analysis for each input seqs.
|
|--Seq01/
|    |-- Seq.fa, Seq.aln, Seq.bed, Seq.dnd
|
|--Seq02/
|    |-- ....
|
|--.../

=head1 AUTHOR

Wang Ming, wangmcas@gmail.com

=cut
