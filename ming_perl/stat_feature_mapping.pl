#!/usr/bin/perl -w

###########################################################################
# Count reads mapping on genomic features: AS,IGR,mRNA,rRNA,tRNA
#
# Using 'featureCounts' to count reads number.
#
# Wang Ming wangmcas(AT)gmail.com 2015-09-09
###########################################################################

use strict;
use warnings;
use Getopt::Std;
use File::Basename qw(basename dirname);
use File::Spec::Functions qw(catfile catdir);
use File::Path qw(make_path remove_tree);
use Cwd qw(abs_path);

### TOOLs
my $find_feature = "/home/wangming/work/bin/temp/find_genomefeatures.pl";
my $featureCounts = "/home/wangming/work/bin/temp/chk_feature2count.pl";
my $sort2bed      = "/home/wangming/work/bin/temp/sort2bed.pl";
# NEED 'featureCounts' in PATH

stat_feature_mapping();

sub stat_feature_mapping {
    my %opts = ();
    getopts("b:f:o:", \%opts);
    usage() if(@ARGV != 1);
    die("[-f] need reference *.fna \n") if(! defined $opts{f});
    my $outdir = $ARGV[0];
    ### prepare features
    my @tmp = prep_feature($opts{f}, $outdir);
    my $feature_dir   = $tmp[0];
    my @feature_files = @{ $tmp[1] };
    ### prepare BAM
    my @bams = readBAMlist($opts{b}, $outdir);
    ### count each features
    my @stats = ();
    for my $f (@feature_files) {
        ### convert bed to gff
        my $gff = $f;
        $gff =~ s/\.bed$/.gff/;
        system("perl $sort2bed -t bed2gff -f exon -i $f -o $gff");
        my @s = count_feature($gff, \@bams, $outdir);
        push @stats, @s;
    }
    ### all
    my $f_stat = catfile($outdir, 'feature.stats');
    open my $fh_f, "> $f_stat" or die "Cannot write $f_stat, $!\n";
    print $fh_f join("\n", @stats) . "\n";
    close $fh_f;
    print join("\n", @stats) . "\n";
}


sub prep_feature {
    my $ref    = $_[0];
    my $outdir = $_[1];
    ### 
    my $fname = $ref;
    $fname =~ s/\.f(n)?a//;
    my $gff = $fname . '.gff';
    my $ptt = $fname . '.ptt';
    my $rnt = $fname . '.rnt';
    for my $f ($gff, $ptt, $rnt) {
        die("[$f] file not exists\n") if(! -e $f);
    }
    ### feature_dir, chr
    my $feature_dir = catdir($outdir, '1.Features');
    make_path($feature_dir) if( ! -d $feature_dir);
    my $chr = 'chr';
    my %h   = ();
    open my $fh_in, "< $gff" or die "Cannot open $gff, $!\n";
    while(<$fh_in>) {
        next if(/^\#|^\s*$/);
        my $n = (split /\t/, $_)[0];
        $h{$n} ++;
        last if($. > 100);
    }
    close $fh_in;
    ### select the name, present in >80% lines
    for my $i (keys %h) {
        if($h{$i} > 80) {
            $chr = $i;
            last;
        }
    }
    ###
    system "perl $find_feature -a $ref -n $chr $feature_dir";
    my @fs = glob("$feature_dir/*.bed");
    return($feature_dir, \@fs);
}

sub readBAMlist {
    my $bamdir = $_[0];
    my $outdir = $_[1];
    my $bamlist = catfile($outdir, 'bam.list');
    system "find $bamdir -name \"*.f.s.bam\" | sort > $bamlist";
    my @bams = ();
    open my $fh_in, "< $bamlist" or die "Cannot open $bamlist, $!";
    while(<$fh_in>) {
        chomp;
        next if(/(^\s*$)|(^\#)/);
        my $b = (split /\s+/)[0];
        die("[$b] in line:$. of $bamlist, not found\n") if(! -e $b);
        push @bams, abs_path($b);
    }
    close $fh_in;
    return @bams;
}

sub count_feature {
    my $file   = $_[0];
    my @bams   = @{ $_[1] };
    my $outdir = $_[2];
    my $count_dir = catdir($outdir, '2.count');
    my $stat_dir  = catdir($outdir, 'output');
    make_path($count_dir);
    make_path($stat_dir);
    my @stats = ();
    for my $bam (@bams) {
        my $stat_sum = featureCounts_count($bam, $file, $count_dir);
        push @stats, $stat_sum;
    }
    my $stat_file = catfile($stat_dir, basename($file) . '.stat');
#    open my $fh_stat, "> $stat_file" or die "Cannot open $stat_file, $!";
#    print $fh_stat join("\n", @stats) . "\n";
#    close $fh_stat;
    return(@stats);
}

sub featureCounts_count {
    my $bam       = $_[0];
    my $file      = $_[1]; # gff
    my $count_dir = $_[2];
    ###
    my $bam_name = basename($bam);
    $bam_name =~ s/(\.|\.s.|\.f.s.|\.trim\.gz\.f\.s\.)bam//;
    my ($file_name) = basename($file) =~ /(.*)\.(\w+)$/; # gff
    my $fcount = catfile($count_dir, $file_name . '_count.txt');
    my $flog   = catfile($count_dir, 'count.log');
    my $para  = '';
    if($bam_name =~ /\_[12]$/) {
    ### PE: dUTP -s 2 (reverse)
        $para = join(' ', '-M --fraction --donotsort -O -f -T 5 -g gene_id -t exon -p -P -d 40 -D 500 -s 2 -a', $file, '-o', $fcount);
    }else {
        $para = join(' ', '-M --fraction --donotsort -O -f -T 5 -g gene_id -t exon -s 1 -a', $file, '-o', $fcount);
    }
    my $run = join(' ', 'featureCounts', $para, $bam, "> $flog 2>&1");
    system"$run";
    my $summary = $fcount . '.summary';
    open my $fh_sum, "< $summary" or die "Cannot open $summary, $!\n";
    my @nums = ();
    while(<$fh_sum>) {
        chomp;
        my $num = (split /\s+/, $_)[1];
        push @nums, $num;
    }
    close $fh_sum;
    my $stat_sum = join("\t", $file_name, $bam_name, @nums);
    return($stat_sum);
}

sub usage {
    die("
Usage: stat_feature_mapping.pl [options] 

Options: -b  <STR>          : the directory of BAM files
         -f  <STR>          : reference file (*.fna), 
                              need *.gff, *.ptt, *.rnt in the same dir
\n");
}

### END OF FILE ###
