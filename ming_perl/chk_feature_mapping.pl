#!/usr/bin/perl -w

###########################################################################
# Count reads mapping on genomic features: AS,IGR,mRNA,rRNA,tRNA          #
#                                                                         #
# Using 'featureCounts' to count reads number.                            #
#                                                                         #
# Wang Ming wangmcas(AT)gmail.com 2015-09-09                              #
###########################################################################

use strict;
use warnings;
use Getopt::Std;
use Cwd qw(abs_path);
use File::Path qw(make_path remove_tree);
use File::Basename qw(basename dirname);
use File::Spec::Functions qw(catfile catdir);

### Requirements ###
# bedtools, featureCounts in the PATH #

### TOOLs
my $find_feature = "/home/wangming/work/bin/temp/getGenomeFeatures.pl";
my $sort2bed      = "/home/wangming/work/bin/temp/sort2bed.pl";
# NEED 'featureCounts' in PATH

stat_feature_mapping();

sub stat_feature_mapping {
    my %opts = (o => 'stat');
    getopts("f:o:n:", \%opts);
    usage() if(@ARGV != 1); # bam_file list, each file in one line.
    die("[-f] *.fna required\n") if(! defined $opts{f});
    my $fna     = $opts{f};
    my $outdir  = $opts{o};
    my $bamlist = $ARGV[0];
    ### prepare working dirs
    my $feature_dir = catfile($outdir, "1.features");
    my $count_dir   = catfile($outdir, "2.count");
    my $output_dir  = catfile($outdir, "output");
    remove_tree($count_dir) if( -d $count_dir);
    make_path($feature_dir) if( ! -d $feature_dir);
    make_path($count_dir)   if( ! -d $count_dir);
    make_path($output_dir)  if( ! -d $output_dir);
    ### prepare features
    my $chr_name = parse_chrname($fna);
    $chr_name    = $opts{n} if(defined $opts{n}); ### specify the chr_name
    my @feature_beds = @{ prep_feature($fna, $chr_name, $feature_dir) };
    ### prepare BAM
    my @bams = readBAMlist($bamlist);
    ### count each features
    my @stats = ();
    for my $f (@feature_beds) {
        ### convert bed to gff
        my $gff = $f;
        $gff =~ s/\.bed$/.gff/;
        system("perl $sort2bed -t bed2gff -f exon -i $f -o $gff");
        my @s = count_feature($gff, \@bams, $count_dir);
        push @stats, @s;
    }
    ### report
    my $f_stat = catfile($output_dir, 'feature.stats');
    open my $fh_f, "> $f_stat" or die "Cannot write $f_stat, $!\n";
    my @outputs = format_output(\@stats);
    print $fh_f join("\n", @outputs) . "\n";
    print join("\n", @outputs) . "\n";
    close $fh_f;
}

sub prep_feature {
    my $ref         = $_[0];
    my $chr_name    = $_[1];
    my $feature_dir = $_[2];
    ### 
    my $fname = $ref;
    $fname =~ s/\.f(n)?a//; # need fna/gff/ptt/rnt files in the same folder
    my $gff = $fname . '.gff';
    my $ptt = $fname . '.ptt';
    my $rnt = $fname . '.rnt';
    for my $f ($gff, $ptt, $rnt) {
        die("[$f] file required\n") if(! -e $f);
    }
    $chr_name =~ s/\|/\\\|/g;
    my $f_stat = `perl $find_feature -f $ref -n $chr_name $feature_dir`;
    my @fs = glob("$feature_dir/*.bed");
    return(\@fs);
}

### get the chr name
sub parse_chrname {
    my $fna = $_[0];
    my $chr_name;
    open my $fh_fna, "< $fna" or die "Cannot open fna, $fna, $!\n";
    while(<$fh_fna>) {
        chomp;
        if(/^\>/) {
            $chr_name = (split /\s+/, $_)[0];
            if($chr_name =~ /\|/) {
                $chr_name = (split /\|/, $chr_name)[3];
            }else {
                next;
            }
        }
    }
    close $fh_fna;
    return($chr_name);
}

### specify the bam files, from the arguments
sub readBAMlist {
    my $bamlist = $_[0];
    my @bams    = ();
    open my $fh_b, "< $bamlist" or die "Cannot open file, $bamlist, $!\n";
    while(<$fh_b>) {
        chomp;
        my $b = (split /\,|\s+/, $_)[0];
        die("[$b] bam file not exists.\n") if( ! -e $b );
        push @bams, abs_path($b);
    }
    close $fh_b;
    return @bams;
}

sub count_feature {
    my $file        = $_[0]; #gff
    my @bams        = @{ $_[1] };
    my $count_dir   = $_[2];
    my ($file_name) = basename($file) =~ /(.*)\.(\w+)$/;
    my $f_count     = catfile($count_dir, $file_name . "_count.txt");
    my $f_summary   = catfile($count_dir, $file_name . "_count.summary");
    my $f_log       = catfile($count_dir, "count.log");
    ### report
    my @sub_counts   = ();
    my @sub_summarys = ();
    my @stats = ();
    my $counter = 1;
    for my $bam (@bams) {
        my $bam_mapped  = qx(samtools idxstats $bam | head -n1 |awk '{print \$3}');
        my $ratio       = sprintf"%.4f", 1000000/$bam_mapped;
        my $sub_count   = catfile($count_dir, "tmp." . $file_name . "." . $counter);
        my $sub_summary = $sub_count . ".summary";
        my $sub_stats   = featureCounts_count($bam, $file, $sub_count, $f_log);
        push @stats, $sub_stats;
        $counter ++;
        calTPM($sub_count, $ratio);
        push @sub_counts, $sub_count;
        push @sub_summarys, $sub_summary;
    }
    ### report summary
    my $run_summary = join(" ", "cat", @sub_summarys, ">", $f_summary);
    system("$run_summary"); ## summary
    ### report count
    combine_count(\@sub_counts, $f_count, $file);
    clear_tmp($count_dir);
    return(@stats);
}

sub clear_tmp {
    my $count_dir = $_[0];
    my @tmp_files = glob("$count_dir\/tmp.*");
    for my $t (@tmp_files) {
        unlink($t);
    }
}

sub calTPM {
    my $count = $_[0];
    my $ratio = $_[1];
    my @tpms  = ();
    open my $fh_c, "< $count" or die "Cannot open $count, $!\n";
    while(<$fh_c>) {
        chomp;
        next if(/^\#|^\s*$/);
        next if(/^Geneid\tChr/);
        my $num = (split /\s+/, $_)[-1];
        my $tpm = sprintf"%.2f", $num * $ratio;
        push @tpms, $_ . "\t" . $tpm;
    }
    close $fh_c;
    open my $fh_o, "> $count" or die "Cannot open $count, $!\n";
    print $fh_o join("\n", @tpms) . "\n";
    close $fh_o;
}

sub combine_count {
    my @sub_counts = @{ $_[0] };
    my $f_count    = $_[1];
    my $file       = $_[2];
    my %nums = ();
    for my $c (sort @sub_counts) {
        open my $fh_c, "< $c" or die "Cannot open $c, $!\n";
        while(<$fh_c>) {
            chomp;
            next if(/^\#|^\s*$/);
            next if(/^Geneid\tChr/);
            my ($id, $num, $tpm) = (split /\s+/, $_)[0, -2, -1];
            push @{$nums{$id}}, ($num, $tpm);
        }
        close $fh_c;
    }
    open my $fh_f, "< $file" or die "Cannot open $file, $!\n";
    open my $fh_o, "> $f_count" or die "Cannot open $f_count, $!\n";
    while(<$fh_f>) {
        chomp;
        next if(/^\#|^\s*$/);
        my ($chr, $start, $end, $strand, $note) = (split /\t/, $_)[0, 3, 4, 6, 8];
        my ($id) = $note =~ /gene_id=(.*)\;locus_tag/;
        my $length = $end - $start + 1;
        print $fh_o join("\t", $id, $chr, $length, $start, $end, $strand,  @{$nums{$id}}) . "\n";
    }
    close $fh_f;
    close $fh_o;
}

sub featureCounts_count {
    my $bam       = $_[0];
    my $file      = $_[1];
    my $f_count   = $_[2];
    my $f_log     = $_[3];
    ###
    my $para  = '';
    my $bam_name = basename($bam);
    $bam_name =~ s/(\.|\.s.|\.f.s.|\.trim\.gz\.f\.s\.)bam//;
    if($bam_name =~ /\_[12]$/) {
    ### PE: dUTP -s 2 (reverse)
        $para = join(' ', '-M --fraction --donotsort -O -f -T 5 -g gene_id -t exon -p -P -d 40 -D 500 -s 2 -a', $file, '-o', $f_count);
    }else {
    ### SE:
        $para = join(' ', '-M --fraction --donotsort -O -f -T 5 -g gene_id -t exon -s 1 -a', $file, '-o', $f_count);
    }
    my $run = join(' ', 'featureCounts', $para, $bam, ">> $f_log 2>&1");
    system"$run";
    ### run commands
    my $summary = $f_count . '.summary';
    open my $fh_sum, "< $summary" or die "Cannot open $summary, $!\n";
    ###
    my $assign     = 0;
    my $not_assign = 0;
    my $unmapped   = 0;
    while(<$fh_sum>) {
        chomp;
        my ($type, $num) = (split /\s+/, $_)[0, 1];
        next if($type eq 'Status');
        if($type eq "Assigned") {
            $assign = $num;
            next;
        }elsif($type eq "Unassigned_Unmapped") {
            $unmapped = $num;
            next;
        }else {
            $not_assign += $num;
        }
    }
    close $fh_sum;
    my ($file_name) = basename($file) =~ /(.*)\.(\w+)$/;
    my $stat_sum = join("\t", $file_name, $bam_name, $assign, $not_assign, $unmapped);
    return($stat_sum);
}

### format the ouput data
sub format_output {
    my @stats = @{ $_[0] };
    my %fmt   = ();
    my %ft    = ();
    for my $s (@stats) {
        my ($feature, $bam, $assign, $not_assign, $unmapped) = split /\s+/, $s;
        $ft{$feature} ++;
        $fmt{$bam}->{map}->{$feature} = $assign;
        $fmt{$bam}->{un}->{$feature}  = $not_assign;
    }
    ###
    my @output   = ();
    my @features = sort keys %ft;
    my $header   = join("\t", "bam_file", @features, @features);
    push @output, $header;
    for my $b (sort keys %fmt) {
        my @lines = ($b);
        for my $f1 (sort keys %{$fmt{$b}->{map}}) {
            push @lines, $fmt{$b}->{map}->{$f1};
        }
        for my $f2 (sort keys %{$fmt{$b}->{un}}) {
            push @lines, $fmt{$b}->{un}->{$f2};
        }
        push @output, join("\t", @lines);
#        print join("\t", @lines) . "\n";
    }
    return(@output);
}

sub usage {
    die("
Usage: stat_feature_mapping.pl [options] <bam_list>

Options: -o  <STR>          : the directory of output files
         -f  <STR>          : reference file (*.fna), 
                              need *.gff, *.ptt, *.rnt in the same dir
         -n  <STR>          : the name of the chr

Example:
perl $0 -b alignment -f ref.fna out_dir
\n");
}

__END__

This script will calculate the read number for genomic features.
1. total number of mapped reads is calculated by samtools/featureCounts;
2. Calculate mRNA,tRNA,rRNA,asmRNA,astRNA,asrRNA reads mapping;
3. Calculate igr mapping by: total_mapped - mRNA - tRNA - rRNA - asmRNA - astRNA - asrRNA;

### END OF FILE ###
