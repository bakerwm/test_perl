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

my $ref_index;
my $ref_gff;
my $out_path;
my $wk_ref;
my @reads;
my $ex_para = '';
my $core = 6;

run_reads();

sub parse_para{
    my %opts = ();
    getopts("d:g:p:o:", \%opts);
    die("Usage: $0 [-d] <ref.idx> [-g] <ref.gff> [-o] <out.dir> reads.list\n") if(@ARGV != 1);
    die("[-d] reference fasta not found.\n") unless defined $opts{d};
    die("[-g] annotation GFF not found.\n") unless defined $opts{g};
    die("[-o] need specify output dir.\n") unless defined $opts{o};
#("[-p] extra parameters for tophat.\n");
    $ex_para = $opts{p} if defined $opts{p};
    my $reads_list = shift(@ARGV);
    die("[$reads_list] file not exists\n") unless -e $reads_list;
    open my $fh_fa, $reads_list or die("$!");
    while(<$fh_fa>){
        chomp;
        next if(/^\s*$/);
        push @reads, $_;
    }
    close $fh_fa;
    $ref_index = $opts{d};
    $ref_gff = $opts{g};
    $out_path = $opts{o};
    # tools
    die('[samtools] not found in $PATH.'."\n") unless(which('samtools'));
    die('[tophat] not found in $PATH.'. "\n") unless(which('tophat'));
    my @index_tmp = glob"$opts{d}*.bt2";
    die('[-d] not found build-index'."\n") if(@index_tmp < 6);
    make_path($out_path) unless -d $out_path;
}

sub run_reads{
    parse_para();
    my @sys_runs = ();
    foreach my $r (sort @reads){
        my $fm = (split /\,|\s+/, $r)[0]; # parse the file name of read1
        my $r_name = basename($fm);
        $r_name =~ s/\.f[astq]+//;
        my $r_outdir = catdir($out_path, $r_name);
        my @runs = perform_tophat($r, $r_outdir);
        push @sys_runs, @runs;
    }
    foreach my $run (@sys_runs){
        print $run, "\n";
        system "$run";
    }
}

sub perform_tophat{
    my ($in, $dir) = @_;
    my $dir_name = basename($dir);
    my $bam_out = catfile($dir, "accepted_hits.bam");
    my $bam_filt = catfile($dir, "$dir_name\.f.bam");
    my $bam_sorted = catfile($dir, "$dir_name\.f.s");
    my @sys_runs = ();
    push @sys_runs, join(' ', ('tophat -p', $core, $ex_para, '-G', $ref_gff, '-o', $dir, '--no-novel-juncs', $ref_index, $in ));
    push @sys_runs, join (' ', ('samtools view -b -F 256', $bam_out, '-o', $bam_filt));
    push @sys_runs, join(' ', ('samtools sort', $bam_filt, $bam_sorted));
    push @sys_runs, join(' ', ('samtools index'," $bam_sorted\.bam"));
    return @sys_runs;    
}

