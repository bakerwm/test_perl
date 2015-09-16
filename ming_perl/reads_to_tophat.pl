#!/usr/bin/perl -w

##############################
# Mapping reads to reference
# genome using tophat
#
# Wang Ming 
# 2015-05-06
##############################

use strict;
use warnings;
use File::Path qw(make_path remove_tree);
use File::Basename qw(basename dirname);
use File::Spec::Functions qw(catfile catdir);
use File::Which;
use Getopt::Std;

my $ex_para = '';
my $core = 12;

reads_to_tophat();
exit(1);

###
sub reads_to_tophat {
    my %opts = ();
    getopts("d:g:p:o:", \%opts);
    usage() if(@ARGV != 1);
    die("[-d] reference fasta not found.\n") unless defined $opts{d};
    die("[-g] annotation GFF not found.\n") unless defined $opts{g};
    die("[-o] need specify output dir.\n") unless defined $opts{o};
    $ex_para = $opts{p} if defined $opts{p};
    my $reads_list = shift(@ARGV);
    my $ref_index  = $opts{d};
    my $ref_gff    = $opts{g};
    my $outdir     = $opts{o};
    make_path($outdir) if(! -d $outdir);
    check_tools();
    check_index($ref_index);
    #
    my @reads = parse_reads($reads_list);
    my @runs  = ();
    for my $read (@reads) {
        my @lines  = get_tophat_cmd($ref_index, $ref_gff, $outdir, $read);
        push @runs, @lines;
    }
    for my $run (@runs) {
        print $run . "\n";
        system "$run";
    }
}

sub parse_reads {
    my $in = shift(@_);
    my @reads = ();
    open my $fh_in, "< $in" or die "Cannot open $in, $!\n";
    while(<$fh_in>) {
        chomp;
        next if(/^\s*$|^\#/);
        push @reads, $_;
    }
    close $fh_in;
    return @reads;
}

sub check_index { # check the bowtie2 index
    my $in = shift(@_);
    my @index_list = glob("$in\*.bt2");
    die('[-d] build-index not found'."\n") if(@index_list < 6);
}

sub get_tophat_cmd {
    my ($index, $gff, $outdir, $read) = @_;
    my $sub_read  = (split /\,|\s+/, $read)[0];
    my ($read_name) = basename($sub_read) =~ /(.*)\.f[astq]+/;
    $read_name   =~ s/\.trim// if($read_name =~ /trim/); # trim the file name
    my $sub_dir  = catdir($outdir, $read_name);
    my $bam_out  = catfile($sub_dir, 'accepted_hits.bam');
    my $bam_filt = catfile($sub_dir, $read_name . '.f.bam');
    my $bam_sort = catfile($sub_dir, $read_name . '.f.s');
    my @runs = ();
    my $libtype = '';
    $libtype = '--library-type fr-secondstrand' if($read =~ /\w\s+\w/); # multiple reads in one line
    push @runs, join(" ", "tophat -p", $core, $libtype, "-G", $gff, "-o", $sub_dir, "--no-novel-juncs", $index, $read);
    push @runs, join(' ', ('samtools view -b -F 256', $bam_out, '-o', $bam_filt));
    push @runs, join(' ', ('samtools sort', $bam_filt, $bam_sort));
    push @runs, join(' ', ('samtools index', " $bam_sort\.bam"));
    return @runs;
}

sub check_tools { # samtools, tophat
    die('[samtools] not found in $PATH.'."\n") unless(which('samtools'));
    die('[tophat] not found in $PATH.'. "\n") unless(which('tophat'));
}

sub usage {
    die("Usage: reads_to_tophat.pl [Options] <in.list>

Options: -d <STR>   : dir of the bowtie2-index 
         -g <STR>   : path to GFF/GTF file
         -o <STR>   : Output dir
         -p <STR>   : Possible additional parameters to tophat.
                      '--library-type fr-secondstrand'

Example:
read_to_tophat.pl -d Ra -g Ra.gff -o Ra_out reads.list
\n");
}

### END OF FILE ###
