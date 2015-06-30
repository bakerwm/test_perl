#!/usr/bin/env perl

##################################
# count reads on each input seq
# using HTSeq-count
##################################

use strict;
use warnings;
use Cwd qw(abs_path cwd);
use File::Which;
use File::Path qw(make_path remove_tree);
use File::Basename qw(dirname basename);
use File::Spec::Functions qw(catdir catfile);
use Getopt::Std;
use Data::Dumper;

use FindBin qw($Bin);
use lib catdir(dirname($Bin), 'lib');

my %func = ('samtools'    => '',
            'bedtools'    => '',
            'htseq-count' => '',
            'sort2bed.pl' => '',
            'featureCounts' => '');

my $f = check_tools(\%func);
%func = %{$f};

&seq2count();

exit(1);

sub seq2count {
    my %opts = (o => 'Count_seqs', p => 'featureCounts', s => 'yes');
    getopts('o:l:f:s:p:', \%opts);

    &usage if(@ARGV == 0);
    die("[-l] Need bam list\n") if(! defined($opts{l}));
    die("[-f] Need reference *.fa n") if(! defined($opts{f}));
    die("[-p: $opts{p}] can be: htseq-count or bedtools.") if(! $opts{p} =~ /^(htseq-count)|(bedtools)|(featureCounts)$/);
    die("[-s: $opts{s}] can be: yes, no or reverse") if(! $opts{s} =~ /^(yes)|(reverse)$/);
    $opts{o} = abs_path($opts{o});
    $opts{f} = abs_path($opts{f});
    make_path($opts{o}) if(! -d $opts{o});

    my $infile = abs_path( shift(@ARGV) );
    my @bams   = readBAMlist($opts{l});
    my $ref    = $opts{f};

    my @runs   = ();
    my $num    = 1;
    my $out_count = catfile($opts{o}, basename($infile));
    $out_count    =~ s/\.txt/_count.txt/;

    # create temp dir
    my $origin_path = abs_path(cwd);
    my $temp_dir = catdir(cwd, 'temp'.int(rand(10000)));
    make_path($temp_dir) if(! -d $temp_dir);
    chdir $temp_dir;
    for my $b (sort @bams) { my $flag = sprintf"%0d", $num; if($opts{p} eq 'htseq-count' ) {
            push @runs, htseq_count($b, $infile, $flag, $opts{s}, $temp_dir);
        }elsif($opts{p} eq 'bedtools') {
            push @runs, bedtools_count($b, $infile, $flag, $opts{s}, $temp_dir);
        }elsif($opts{p} eq 'featureCounts') {
            push @runs, featureCounts_count($b, $infile, $flag, $opts{s}, $temp_dir);
        }
        $num ++;
    }
   
    if($opts{p} eq 'featureCounts') {
        my $tmp_summary = catfile($temp_dir, basename($infile) . '.tmp.summary.*');
        my $out_summary = $out_count;
        $out_summary    =~ s/\.txt/.summary/;
        push @runs, "cat $tmp_summary > $out_summary";
    }
    
    my $tmp2_file = catfile($temp_dir, basename($infile).'.tmp2');
    my $TPM_files = catfile($temp_dir, basename($infile).'.TPM.*');
    push @runs, "sort -k1 $infile -o $infile";
    push @runs, "paste $infile $TPM_files > $out_count";

    for my $r (@runs) {
        print $r . "\n";
        system"$r";
    }
    chdir $origin_path;
    remove_tree($temp_dir);
}

#sub count_seq {
sub htseq_count {
    my ($bam, $in, $num, $strand, $outdir) = @_;
    # determine strand
    my $bam_name = basename($bam);
    $bam_name    =~ s/(\.|\.s.|\.f.s.|\.trim\.gz\.f\.s\.)bam//;
    my $lib_type = '';
    if($strand eq 'yes') {
        $lib_type = ($bam_name =~ /\_[12]$/)?'reverse':'yes'; # PE=reverse, SE=yes
    }elsif($strand eq 'reverse') {
        $lib_type = ($bam_name =~ /\_[12]$/)?'yes':'reverse';
    }else {
        #
    }
    # run 
    my @runs      = ();
    chomp(my $chr = qx($func{'samtools'} view $bam | head -n1 | awk '{print \$3}'));
    my $in_gff    = basename($in);
    $in_gff       =~ s/\.txt$/.gff/;
    $in_gff       = catfile($outdir, $in_gff);
    my $in_tmp    = catfile($outdir, basename($in).'.tmp');
    my $in_tmp2   = catfile($outdir, basename($in).'.tmp2');
    my $in_TPM    = catfile($outdir, basename($in).'.TPM.'.$num);
    if( not_blank_file($in) ) {
        push @runs, "perl $func{'sort2bed.pl'} -t sort2gff -f exon -s $chr -i $in -o $in_gff";
        push @runs, "$func{'htseq-count'} -q -f bam -s $lib_type -t exon $bam $in_gff > $in_tmp";
        push @runs, "sort -k1 $in_tmp | sed -e \'/^\_/d\' > $in_tmp2";
        # count tpm
        my $bam_mapped = qx($func{'samtools'} idxstats $bam | head -n1 |awk '{print \$3}');
        my $ratio = sprintf"%.4f", 1000000/$bam_mapped;
        push @runs, "awk \'{printf(\"\%s\\t\%s\\t\%.4f\\n\", \$1, \$2, \$2*$ratio)}\'  \< $in_tmp2 | cut -f2-3 > $in_TPM";
    }
    return @runs;
}

sub featureCounts_count {
    my ($bam, $in, $num, $strand, $outdir) = @_;
    # determine strand
    my $bam_name = basename($bam);
    $bam_name    =~ s/(\.|\.s.|\.f.s.|\.trim\.gz\.f\.s\.)bam//;
    my $lib_type = '';
    if($strand eq 'yes') {
        $lib_type = ($bam_name =~ /\_[12]$/)?'2':'1'; # PE=2, SE=1
    }elsif($strand eq 'reverse') {
        $lib_type = ($bam_name =~ /\_[12]$/)?'1':'2'; # reverse
    }else {
        #
    }
    my @runs      = ();
    chomp(my $chr = qx($func{'samtools'} view $bam | head -n1 | awk '{print \$3}'));
    my $in_gff    = basename($in);
    $in_gff       =~ s/\.txt$/.gff/;
    $in_gff       = catfile($outdir, $in_gff);
    my $in_tmp    = catfile($outdir, basename($in).'.tmp');
    my $in_tmp2   = catfile($outdir, basename($in).'.tmp2');
    my $in_TPM    = catfile($outdir, basename($in).'.TPM.'.$num);
    if( not_blank_file($in) ) {
        push @runs, "perl $func{'sort2bed.pl'} -t sort2gff -f exon -s $chr -i $in -o $in_gff";
        push @runs, "$func{'featureCounts'} -p -T 5 -s $lib_type -M -O --donotsort -t exon -g gene_id -a $in_gff -o $in_tmp $bam >>$in\.log 2>&1"; # need: -O 
# keep summary
        push @runs, "mv $in_tmp\.summary $in_tmp\.summary\.$num";
        push @runs, "cat $in_tmp | sed  \'1,2 d\' | sort -k1 > $in_tmp2"; 
        # count tpm
        my $bam_mapped = qx($func{'samtools'} idxstats $bam | head -n1 |awk '{print \$3}');
        my $ratio = sprintf"%.4f", 1000000/$bam_mapped;
        push @runs, "awk \'{printf(\"\%s\\t\%s\\t\%.4f\\n\", \$1, \$7, \$7*$ratio)}\'  \< $in_tmp2 | cut -f2-3 > $in_TPM";
    }
    return @runs;
}

sub bedtools_count {
    my ($bam, $in, $num, $strand, $outdir) = @_;
    my $bam_name = basename($bam);
    $bam_name    =~ s/(\.|\.s.|\.f.s.|\.trim\.gz\.f\.s\.)bam//;
    my @runs     = ();
    my $str_type = '';
    if($strand eq 'yes') {
        $str_type = '-s';
    }elsif($strand eq 'reverse') {
        $str_type = '-S';
    }else {
        $str_type = '';
    }
    my $chr    = qx($func{'samtools'} view $bam | head -n1 | awk '{print \$3}');
    chomp($chr);
    my $in_bed = basename($in);
    $in_bed    =~ s/\.txt/.bed/;
    $in_bed    = catfile($outdir, $in_bed);
    my $in_tmp = catfile($outdir, basename($in).'.tmp');
    my $in_TPM = catfile($outdir, basename($in).'.TPM.'.$num);
    if( not_blank_file($in) ) {
        push @runs, "perl $func{'sort2bed.pl'} -t sort2bed -s $chr -i $in -o $in_bed";

#        push @runs, "$func{'bedtools'} intersect -c $str_type -a $in_bed -b $bam  > $in_tmp";
        push @runs, "$func{'bedtools'} multicov -split -D $str_type -bams $bam -bed $in_bed > $in_tmp";
        # count tpm
        my $bam_mapped = qx($func{'samtools'} idxstats $bam | head -n1 |awk '{print \$3}');
        my $ratio      = sprintf"%.4f", 1000000/$bam_mapped;
        push @runs, "awk \'{printf(\"\%s\\t\%s\\t\%.4f\\n\", \$1, \$7, \$7*$ratio)}\'  \< $in_tmp | cut -f2-3 > $in_TPM";
    }
    return @runs;
}

sub readBAMlist {
    my $in = shift(@_);
    die("[$in] BAM list file not found.\n") if(! -e $in);
    my @bams = ();
    open my $fh_in, "< $in" or die "$!";
    while(<$fh_in>) {
        chomp;
        next if(/(^\s*$)|(^\#)/);
        my $b = (split /\s+/)[0];
        die("[$b] in line:$. of $in, not found\n") if(! -e $b);
        push @bams, abs_path($b);
    }
    close $fh_in;
    return @bams;
}

sub not_blank_file {
    my $in = shift(@_);
    if( -e $in) {
        open my $fh_in, "< $in" or die "$!";
        my $num = 0;
        while(<$fh_in>) {
            chomp;
            next if(/(^\s*$)|(^\#)/);
            $num ++;
        }
        close $fh_in;
        return $num;    
    }else {
        return 0;
    }
}

# htseq-count
# samtools
sub check_tools {
    my $tool = shift(@_);
    my %func = %{$tool};
    my @missing = ();
    for my $f (sort keys %func) {
        if( tool_path($f) ) {
            $func{$f} = tool_path($f);
        }else {
            push @missing, $f;
        }
    }
    if(@missing) {
        print 'The following tool(s) missing:' . "\n\n";
        my $err_log = '';
        for my $m (sort @missing) {
            $err_log .= sprintf("%-15s : Not found in \$PATH, \$Bin, or ~/work/bin/temp/\n", $m);
        }
        die("$err_log\n");
    }
    my $flag = (@missing)?0:1;
    return($flag, \%func);
    #
    sub tool_path {
        my $t = shift(@_);
        my $perl_path = $ENV{HOME}. '/work/bin/temp';
        if( -e catfile($perl_path, $t)) {
            return catfile($perl_path, $t);
        }elsif( which($t) ) {
            return $t;
        }elsif( -e catfile($Bin, $t) ) {
            return catfile($Bin, $t);
        }else {
            return 0;
        }
    }
}

sub usage{
    die(qq/
Usage: chk_seq2count.pl [options] <in.txt>

Options: -o     Output path. [Count_seqs]
         -l     The bam list files, only one bam on each line
         -f     The reference Fasta file
         -s     Strand : yes no reverse [yes]
         -p     The program to count reads: bedtools, htseq-count, featureCounts
                [featureCounts]

\n/);
}
