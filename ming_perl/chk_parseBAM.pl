#!/usr/bin/env perl

#######################################################
# parsing the BAM output from bowtie/bowtie2/tophat,
# 1. stat:  stat mapping reads for each sample
# 2. view:  create bedgraph views for each bams
# 3. tags:  find tags in each bam
#
#######################################################

use strict;
use warnings;
use Cwd qw(abs_path cwd);
use File::Which;
use File::Path qw(make_path remove_tree);
use File::Basename qw(basename dirname);
use File::Spec::Functions qw(catfile catdir);
use POSIX qw(strftime);
use Getopt::Std;
use Data::Dumper;

my %func = ("samtools"     => '',
            "bedtools"     => '',
            "htseq-count"  => '',
            "sort2bed.pl"  => '',
            "search_cov_regions.pl"   => '',
            "sort_to_position_sig.pl" => '',
            "sort2candi_v1.pl"        => '',
            "chk_seq2rnaz.pl"         => '',
            "featureCounts"           => '');

&usage if(@ARGV < 1);

my $command = shift(@ARGV);
my %prog = (stat => \&statBAM,
            view => \&viewBAM,
            tags => \&BAM2tags);

die("Unknown command [$command] \n") if (!defined($prog{$command}));

my $fv = check_tools(\%func);
%func  = %{$fv};

&{$prog{$command}};


exit(0);

sub statBAM {
    my %opts = (o => 'BAM_output');
    getopts('o:', \%opts);
    die(qq/
Usage: parseBAM.pl stat [options] <inbam.list>

<inbam.list> each line contain one BAM file

Options: -o    output dir, [BAM_output]
\n/) if (@ARGV == 0);

    my $stat_dir  = catdir($opts{o}, 'stat_bams');
    my $stat_file = catfile($stat_dir, 'bam.stat');
    make_path($stat_dir) if(! -d $stat_dir);
    open my $fh_st, "> $stat_file" or die "Cannot open file: $stat_file, $!";
    my @stats = ();
    for my $bam (sort &readBAMlist) {
        my $mapped = qx($func{'samtools'} idxstats $bam | head -n1 |awk '{print \$3}');
        chomp($mapped);
        push @stats, $mapped;
        print $bam . "\n";
        print $fh_st $bam . ":\t" . $mapped, "\n";
    }
    close $fh_st;
    return @stats;
}

sub viewBAM {
    my %opts = (o => 'BAM_output', s => 1);
    getopts('s:f:o:', \%opts);
    die(qq/
Usage: parseBAM.pl view [options] <inbam.list>

<inbam.list> each line contain one BAM file

Options: -o str     output dir, [BAM_output]
         -s str     convert the bedgraph by scale number, [1]
                    0=(to 1 M total reads), [0-1]=(the input scale)
         -f str     reference file [fa]
Example:

parseBAM.pl view -o outdir -s 1 -f ref.fa inbam.list > log
\n/) if (@ARGV == 0);
    die("[-f] reference file not exist") if(! -e $opts{f});
    system"$func{'samtools'} faidx $opts{f}";
    my $ref_idx  = $opts{f} . '.fai';
    my $view_dir = catdir($opts{o}, 'view_bedgraph');
    my $view_log = catfile($opts{o}, 'view_bedgraph.log');
    open my $fh_view, "> $view_log" or die "$!";
    for my $bam (sort &readBAMlist) {
       my @runs = bam2bg($bam, $ref_idx, $opts{s}, $view_dir);
       for my $r (@runs) {
            print $r . "\n";
            print $fh_view $r . "\n";
            system"$r";
        }
    }
    close $fh_view;
}

sub BAM2tags {
    my %opts = (o => 'BAM_output', 
                s => 1, 
                c => 100, 
                m => 1,
                e => 1, 
                z => 1,
                d => '/home/wangming/work/database/H37Rv/SixRv.fa');
    getopts("f:g:s:c:m:e:z:o:d:", \%opts);
    die(qq/
Usage: parseBAM.pl tags [options] <inbam.list>

<ibam.list> each line contain one BAM file

Options: -o str     output dir, [BAM_output]
         -f str     reference file [fa]
         -g str     annotation file [gff]
         -s float   scale for calling coverage [1]
                    covert graph multiple by scale:
                    0=(to 1 M reads), [0-1]=(the input scale)
         -c INT     cut-off for determine edges of tags [100]
         -m INT     merge tags from all samples, 0=no, 1=yes, [1]
         -e INT     count reads on each tag by htseq-count. very slow. 0=no, 1=yes, [1]
         -z INT     apply RNAz analysis, 0=no, 1=yes, [1]
                    Only report one region with the highest z-score if
                    the input sequence contain multiple regions with 
                    z-score > 0.5. (criteria)
         -d str    The database for RNAz analysis. [SixRv]
                    see: ~\/work\/database\/H37Rv\/SixRv.fa

Example:

parseBAM.pl tags -o outdir -f ref.fa -g ref.gff inbam.list > log
\n/) if(@ARGV == 0);
    die("[-f: $opts{f}] reference file not exist\n") if(! -e $opts{f});
    die("[-g: $opts{g}] annotation file not exist\n") if(! -e $opts{g});
    my $RNAz_db = abs_path($opts{d});
    system"samtools faidx $opts{f}";
    my $ref_idx = $opts{f} . '.fai';
    my $tag_dir = catdir($opts{o}, 'find_tags');
    my @tags_files = ();
# find tags
    for my $bam (sort &readBAMlist) {
        $opts{s} = 0 if(! defined $opts{s});
        my ($tag_out, @run_tags) = bam2tags($bam, $ref_idx, $opts{f}, $opts{g}, $opts{s}, $opts{c}, $tag_dir);
        push @tags_files, $tag_out;
        run_cmd('0', @run_tags);
    }
# merge tags
    my $merge_dir  = catdir($opts{o}, 'merge_tags');
    make_path($merge_dir) if(! -d $merge_dir);
    my @run_merges = merge2tags($merge_dir, @tags_files);
    run_cmd('1', @run_merges);
# filter tags
    my $merged_file = catfile($merge_dir, 'merged.bed');
    my $lib_num     = @tags_files;
#    my $lib_num     = 1;
    my @sub_beds    = ();
    if( not_blank_file($merged_file) ) {
        @sub_beds   = filtermerged($merge_dir, $merged_file, $lib_num);
    }
    for my $tag (sort @sub_beds) {
        my @run_filts  = ();
        my $tag_dir   = dirname($tag);
        my $tag_new   = catfile($tag_dir, 'tag.newID.bed');
        my $tag_txt   = catfile($tag_dir, 'tag.newID.txt');
        my $tag_pos   = catfile($tag_dir, 'tag.newID.pos.txt');
        my $tag_sRNA  = catfile($tag_dir, 'tag.newID.pos_sRNA.txt');
        my $tag_count = catfile($tag_dir, 'tag.newID.pos_sRNA.count.txt');
        renameID($tag, 4, $tag_new); # id in col-4
        push @run_filts, "perl $func{'sort2bed.pl'} -t bed2sort -i $tag_new -o $tag_txt";
        push @run_filts, "perl $func{'sort_to_position_sig.pl'} -f $opts{f} -g $opts{g} $tag_txt > $tag_pos";
        push @run_filts, "perl $func{'sort2candi_v1.pl'} $tag_pos";
        run_cmd('1', @run_filts);
# count + tpm
        my @run_counts = ();
        my @hts = ();
        if($opts{e}) {
            my $flag = 1;
### HTSeq-count
#            my @hts  = ();
#            for my $bam (sort &readBAMlist) {
#                @hts = txt2count($bam, $tag_sRNA, $tag_count, $flag);
#                push @run_counts, @hts;
#                $flag ++;
#            }
#            if(@hts) {
#                push @run_counts, "paste $tag_sRNA $tag_sRNA\.TPM\.* > $tag_count";
#            }
### FeatureCounts
            my @feas = ();
            for my $bam (sort &readBAMlist) {
                 @feas = featureCounts2count($bam, $tag_sRNA, $tag_count, $flag);
                push @run_counts, @feas;
                $flag ++;
            }
            if(@feas) {
                push @run_counts, "paste $tag_sRNA $tag_sRNA\.TPM\.* > $tag_count";
            }
        }
#        push @run_counts, "rm -rf $tag_sRNA\.tmp $tag_sRNA\.tmp2 $tag_sRNA\.TPM\.\* ";
        run_cmd('1', @run_counts);        
# rnaz analysis
        my @run_rnazs = ();
        if($opts{z}) {
            my $rnaz_dir = catdir(dirname($tag_sRNA), 'RNAz_out');
            make_path($rnaz_dir) if (! -d $rnaz_dir);
            push @run_rnazs, seq2RNAz($tag_sRNA, $opts{f}, $RNAz_db, $rnaz_dir);
        }
        run_cmd('1', @run_rnazs);
# wrap output
        my $wrap_dir = catdir($opts{o}, 'report');
        make_path($wrap_dir) if(! -d $wrap_dir);
        my $sRNA_RNAz   = catfile(catdir($tag_dir, 'RNAz_out'), 'best_RNAz.bed');
        my $sRNA_report = catfile($wrap_dir, basename($tag_dir) . '.report.txt');
#exit (1);        
        my $rpt_lines   = wrap_output($tag_txt, $tag_sRNA, $tag_count, $sRNA_RNAz);
        my $header = '#colum name:[1-12]ID,chr,length,start,end,strand,pre-gene,gap-1,next-gene,gap-2,direction,description'.
                     "\n" . '#exp [13...] count:tpm' .'[last 2-col]RNAz old_ID';
        open my $fh_rpt, "> $sRNA_report" or die "$!";
        print $fh_rpt $header, "\n";
        print $fh_rpt $rpt_lines;
        close $fh_rpt;
    }
}

#
sub run_cmd {
    my @in = @_;
    my $run = shift(@in);
    print join("\n", @in). "\n\n";
    if($run) {
        for my $r (@in) {
            system"$r";
        }
    }
}

sub readBAMlist {
    my $bamlist = $ARGV[0];
    die("[$bamlist] file not exists\n") if(! -e $bamlist);
    my @bamlists = ();
    open my $fh_bam, "< $bamlist" or die "Cannot open file $bamlist, $!\n";
    while(<$fh_bam>) {
        next if(/(^\s*$)|^\#/); # blank lines and #comment lines
        s/^\s+|\s+$//;
        die("[$_] file not exists in line-[$.] of $bamlist\n") if(! -e $_);
        push @bamlists, $_;
    }
    close $fh_bam;
    return @bamlists;
}

sub splitPEBAM {
    my ($bam, $fwdbam, $revbam) = @_;
    my @ps  = ();
    push @ps, "$func{'samtools'} view -b -f 128 -F 16 $bam > fwd1.bam";
    push @ps, "$func{'samtools'} view -b -f 80 $bam > fwd2.bam";
    push @ps, "$func{'samtools'} index fwd1.bam";
    push @ps, "$func{'samtools'} index fwd2.bam";
    push @ps, "$func{'samtools'} merge -f fwd.bam fwd1.bam fwd2.bam";
    push @ps, "$func{'samtools'} view -b -f 144  $bam > rev1.bam";
    push @ps, "$func{'samtools'} view -b -f 64 -F 16 $bam > rev2.bam";
    push @ps, "$func{'samtools'} index rev1.bam";
    push @ps, "$func{'samtools'} index rev2.bam";
    push @ps, "$func{'samtools'} merge -f rev.bam rev1.bam rev2.bam";
    push @ps, "mv -f fwd.bam $fwdbam";
    push @ps, "mv -f rev.bam $revbam";
    push @ps, "rm -f fwd1.bam* fwd2.bam* rev1.bam* rev2.bam*";
    return @ps;
}

sub splitSEBAM {
    my ($bam, $fwdbam, $revbam) = @_;
    my @ps = ();
    push @ps, "$func{'samtools'} view -b -F 16 $bam > $fwdbam";
    push @ps, "$func{'samtools'} view -b -f 16 $bam > $revbam";
    return @ps;
}

sub bam2bg {
    my ($bam, $ref_idx, $s, $outdir) = @_;
    my $t_mapped = qx($func{'samtools'} idxstats $bam | head -n1 |awk '{print \$3}');
    chomp($t_mapped);
    die("It's a blank BAM file") if(! $t_mapped);
    my $scale = sprintf"%.4f", 1000000/$t_mapped;
    $scale = $s if ($s > 0);
    my $bam_name = basename($bam);
    $bam_name =~ s/(\.|\.s.|\.f.s.|\.trim\.gz\.f\.s\.)bam//;
    my $smp_dir = catdir($outdir, $bam_name);
    make_path($smp_dir) if(! -d $smp_dir);
    my $fwd_bam = catfile($smp_dir, $bam_name.'.fwd.bam');
    my $rev_bam = catfile($smp_dir, $bam_name.'.rev.bam');
    my @runs = ();
    if($bam_name =~ /_[12]$/) {
        @runs = splitPEBAM($bam, $fwd_bam, $rev_bam);
    }else{
        @runs = splitSEBAM($bam, $fwd_bam, $rev_bam);
    }
    my $fwd_bg = catfile($smp_dir, $bam_name . '.fwd.bedgraph');
    my $rev_bg = catfile($smp_dir, $bam_name . '.rev.bedgraph');
    push @runs, "bedtools genomecov -bg -split -scale $scale -ibam $fwd_bam -g $ref_idx > $fwd_bg";
    push @runs, "bedtools genomecov -bg -split -scale $scale -ibam $rev_bam -g $ref_idx > $rev_bg";
    return @runs;
}

sub bam2tags {
    my ($bam, $ref_idx, $ref, $gff, $s, $cov_cutoff, $outdir) = @_;
    my $t_mapped = qx($func{'samtools'} idxstats $bam | head -n1 |awk '{print \$3}');
    chomp($t_mapped);
    die("It's a blank BAM file") if(! $t_mapped);
    my $scale = sprintf"%.4f", 1000000/$t_mapped;
    $scale = $s if ($s > 0);
    my $bam_name = basename($bam);
    $bam_name =~ s/(\.|\.s.|\.f.s.|\.trim\.gz\.f\.s\.)bam//;
    my $smp_dir = catdir($outdir, $bam_name);
    make_path($smp_dir) if(! -d $smp_dir);
    my $fwd_bam = catfile($smp_dir, $bam_name.'.fwd.bam');
    my $rev_bam = catfile($smp_dir, $bam_name.'.rev.bam');
    my @runs = ();
    if($bam_name =~ /_[12]$/) {
        @runs = splitPEBAM($bam, $fwd_bam, $rev_bam);
    }else{
        @runs = splitSEBAM($bam, $fwd_bam, $rev_bam);
    }
    my $fwd_cov   = catfile($smp_dir, $bam_name . '.coverage.p');
    my $rev_cov   = catfile($smp_dir, $bam_name . '.coverage.n');
    my $tag_p     = catfile($smp_dir, $bam_name . '.tag.p');
    my $tag_n     = catfile($smp_dir, $bam_name . '.tag.n');
    my $tag       = catfile($smp_dir, $bam_name . '.tag.txt');
    my $tag_pos   = catfile($smp_dir, $bam_name . '.tag.pos.txt');
    push @runs, "$func{'bedtools'} genomecov -d -split -scale $scale -ibam $fwd_bam -g $ref_idx > $fwd_cov";
    push @runs, "$func{'bedtools'} genomecov -d -split -scale $scale -ibam $rev_bam -g $ref_idx > $rev_cov";
    push @runs, "perl $func{'search_cov_regions.pl'} -c $cov_cutoff -s + $fwd_cov > $tag_p";
    push @runs, "perl $func{'search_cov_regions.pl'} -c $cov_cutoff -s - $rev_cov > $tag_n";
    push @runs, "cat $tag_p $tag_n > $tag";
    push @runs, "perl $func{'sort_to_position_sig.pl'} -f $ref -g $gff $tag > $tag_pos";
    push @runs, "perl $func{'sort2candi_v1.pl'} $tag_pos";
    return ($tag, @runs);
}

sub merge2tags {
    my $outdir  = shift(@_);
    my @infiles = @_;
    my @runs = ();
    my $count = 1;
    my $beds_line = '';
    for my $i (sort @infiles) {
        my $i_bed = catfile($outdir, basename($i));
        $i_bed    =~ s/\.txt$/.bed/;
        # add prefix to the id
        my $flag = sprintf"%02d", $count;
        if(@infiles == 1) {
            push @runs, "sed -e \'s/^/LibN\_/\' $i > $i\.tmp";
        }else{
            push @runs, "sed -e \'s/^/Lib$flag\_/\' $i > $i\.tmp";
        }
        # sort 2 bed
        push @runs, "perl $func{'sort2bed.pl'} -t sort2bed -i $i\.tmp -o $i_bed";

        push @runs, "rm -rf $i\.tmp";
        $beds_line .= $i_bed. " ";
        $count++;
    }
    my $bed_all   = catfile($outdir, 'all_tags.bed');
    my $bed_merge = catfile($outdir, 'merged.bed');
    push @runs, "cat $beds_line | sort -k1,1 -k2,2n | cut -f1-6 > $bed_all";
    push @runs, "$func{'bedtools'} merge -s -d -1 -c 4,5,6 -o distinct,distinct,distinct -i $bed_all > $bed_merge";
    push @runs, "rm -rf $beds_line";
    return @runs;
}

sub renameID {
    my ($in, $col, $out) = @_;
    $col --; # perl is 0-leftmost index
    my $newline = '';
    if( not_blank_file($in) ) {
        open my $fh_in, "< $in" or die "$!";
        while(<$fh_in>) {
            chomp;
            next if(/(^\s*$)|(^\#)/);
            my @tabs    = split /\t/;
            my $newid   = my $id = $tabs[$col];
            $newid      = (split /\,|\:/, $newid)[0];
            $tabs[$col] = $newid;
            $newline   .= join("\t",@tabs, $id) . "\n";
        }
        close $fh_in;
    }
    open my $fh_out, "> $out" or die "$!";
    print $fh_out $newline;
    close $fh_out;
}

sub filtermerged {
    # it's designed for tags merged from differert libraries in various length.
    # default: lib01 18-40 nt, lib02 40-80 nt, lib03 80-140 nt, lib04 >140 nt
    # with -40 to +40 range of lib
    my ($outdir, $merged, $files_num) = @_;
    my %lib = ();
    my %tag = ();
    my $count = 1;
    # create sub dirs
    for(my $i=1; $i<=$files_num; $i++){
        my $flag = sprintf "%02d", $i;
        my $sub_dir = catdir($outdir, "Lib$flag");
        $sub_dir    = catdir($outdir, 'LibN') if($files_num == 1);
        make_path($sub_dir);
        $lib{'tag'}->{$count} = catfile($sub_dir, 'tag.bed');
        $lib{'del'}->{$count} = catfile($sub_dir, 'del.bed');
        $count ++;
    }
    $lib{'others'}->{0}   = catfile($outdir, 'unfiltered.txt');
    @{$tag{'other'}->{0}} = ();
    open my $fh_mg, "< $merged" or die "Cannot open $merged, $!";
    while(<$fh_mg>) {
        chomp;
        next if(/(^\s*$)|(^\#)/);
        my $line = $_;
        my ($start, $end) = (split /\t/, $_)[1,2];
        my $len = $end - $start + 1;
        if(/Lib04/){
            if($len >=100){
                push @{$tag{'tag'}->{4}}, $line;
            }else{
                push @{$tag{'del'}->{4}}, $line;
            }
        }elsif(/Lib03/){
            if($len >= 40 && $len <= 180){
                push @{$tag{'tag'}->{3}}, $line;
            }else{
               push @{$tag{'del'}->{3}}, $line;
            }
        }elsif(/Lib02/){
            if($len >= 40 && $len <= 120){
                push @{$tag{'tag'}->{2}}, $line;
            }else{
                push @{$tag{'del'}->{2}}, $line;
            }
        }elsif(/Lib01/){
            if($len >= 20 && $len <= 80){
                push @{$tag{'tag'}->{1}}, $line;
            }else{
                push @{$tag{'del'}->{1}}, $line;
            }
        }elsif(/LibN/){
            if($len >= 20){
                push @{$tag{'tag'}->{1}}, $line;
            }else{
                push @{$tag{'del'}->{1}}, $line;
            }
        }else{
            push @{$tag{'other'}->{0}}, $line;
        }
    }
    close $fh_mg;
    for my $type (sort keys %lib){
        for my $n (sort keys %{$lib{$type}}){
            open my $fh_n, "> $lib{$type}->{$n}" or die "$!";
            if(exists $tag{$type}->{$n}){
                print $fh_n join("\n", @{$tag{$type}->{$n}}), "\n";
            }else{
            }
            close $fh_n;
        }
    }
    return (sort values %{$lib{'tag'}});
}

sub txt2count {
    my ($bam, $infile, $outfile, $flag) = @_;
    my $bam_name = basename($bam);
    $bam_name   =~ s/(\.|\.s.|\.f.s.|\.trim\.gz\.f\.s\.)bam//;
    my $lib_type = ($bam_name =~ /\_[12]$/)?'reverse':'yes'; # PE=reverse, SE=yes
    my @runs = ();
    my $infile_gff  = $infile;
    $infile_gff =~ s/\.txt$/.gff/;
    if( not_blank_file($infile) ) {
        push @runs, "perl $func{'sort2bed.pl'} -t sort2gff -f exon -i $infile -o $infile_gff";
        push @runs, "$func{'htseq-count'} -q -f bam -s $lib_type -t exon $bam $infile_gff > $infile\.tmp";
        push @runs, "sort -k1 $infile.tmp | sed -e \'/^\_/d\' > $infile\.tmp2";
# add tpm
        # count tpm
        my $mapped  = qx($func{'samtools'} idxstats $bam | head -n1 |awk '{print \$3}');
        my $m_scale = sprintf"%.4f", $mapped/1000000;
        push @runs, "awk \'{printf(\"\%s\\t\%s\\t\%.4f\\n\", \$1, \$2, \$2/$m_scale)}\' $infile\.tmp2 | cut -f2-3 > $infile\.TPM\.$flag";
    }
    return @runs;

}

sub featureCounts2count {
    my ($bam, $infile, $outfile, $flag) = @_;
    my $bam_name = basename($bam);
    $bam_name   =~ s/(\.|\.s.|\.f.s.|\.trim\.gz\.f\.s\.)bam//;
    my $lib_type = ($bam_name =~ /\_[12]$/)?'2':'1'; # PE=2, SE=1
    my @runs = ();
    my $infile_gff = $infile;
    $infile_gff =~ s/\.txt/.gff/;
    if( not_blank_file($infile) ) {
        push @runs, "perl $func{'sort2bed.pl'} -t sort2gff -f exon -i $infile -o $infile_gff";
        push @runs, "$func{'featureCounts'} -p -T 5 -s $lib_type -M --donotsort -t exon -g gene_id -a $infile_gff -o $infile\.count\.$flag $bam >$infile\.log 2>&1";
        push @runs, "sed  \'1,2 d\' $infile\.count\.$flag > $infile\.tmp2\.$flag";
        # count TPM
        my $mapped  = qx($func{'samtools'} idxstats $bam | head -n1 |awk '{print \$3}');
        my $m_scale = sprintf"%.4f", $mapped/1000000;
        push @runs, "awk \'{printf(\"\%s\\t\%s\\t\%.4f\\n\", \$1, \$7, \$7/$m_scale)}\' $infile\.tmp2\.$flag | cut -f2-3 > $infile\.TPM\.$flag";
    }
    return @runs;
}

sub seq2RNAz {
    my ($txt, $ref, $rnaz_db, $outdir) = @_;
    my $txt_fa = $txt;
    $txt_fa    =~ s/\.txt$/.fa/;
    my $rnaz_log = catfile($outdir, 'rnaz.log');
#    my $rnaz_bed = catfile($outdir, 'best_RNAz.bed');
    my @runs = ();
    push @runs, "perl $func{'sort2bed.pl'} -t sort2fa -g $ref -i $txt -o $txt_fa";
    push @runs, "perl $func{'chk_seq2rnaz.pl'} RNAz -d $rnaz_db -o $outdir $txt_fa > $rnaz_log 2>&1 ";
    return @runs;
}

sub wrap_output {
    my ($info, $pos, $exp, $z) = @_;
    my $vf = fetch_id($info);
    my $vz = fetch_id($z);
    my $ve = fetch_count($exp);
    my %hf = %{$vf};
    my %hz = %{$vz};
    my %he = %{$ve};

    my $rpt_out = '';
    if( not_blank_file($pos) ) {
        open my $fh_pos, "< $pos" or die "$!";
        while(<$fh_pos>) {
            chomp;
            next if(/(^\s*$)|(^\#)/);
            my $id = (split /\t/)[0];
            my $zscore = (exists $hz{$id})?$hz{$id}:'-';
            my $note   = (exists $hf{$id})?$hf{$id}:'-';
            my $exp    = (exists $he{$id})?$he{$id}:'';
            $rpt_out  .= join("\t", $_, $exp, $zscore, $note) . "\n";
        }
        close $fh_pos;
        return $rpt_out;
    }else {
        return '';
    }
}

sub fetch_id {
    my $in   = shift(@_);
    my %info = ();
    if( not_blank_file($in) ) {
        open my $fh_in, "< $in" or die "Cannot open $in, $!\n";
        while(<$fh_in>) {
            chomp;
            next if(/(^\s*$)|(^\#)/);
            my ($id, $note) = (split /\t/)[0, -1];
            if($in =~ /RNAz\.bed$/) {
                $id =~ s/^[a-zA-Z0-9]+\_//;
            }
            $info{$id} = $note;
        }
        close $fh_in;
    }
    return \%info;
}

sub fetch_count {
    my $in   = shift(@_);
    my %info = ();
    if( not_blank_file($in) ) {
        open my $fh_in, "< $in" or die "Cannot open $in, $!\n";
        while(<$fh_in>) {
            chomp;
            next if(/(^\s*$)|(^\#)/);
            my @tabs = split /\t/, $_, 13;
            $info{$tabs[0]} = $tabs[-1];
        }
        close $fh_in;
    }
    return \%info;
}

sub not_blank_file {
    my $in = shift(@_);
    if( -e $in ) {
        open my $fh_in, "< $in" or die "$!";
        my $count = 0;
        while(<$fh_in>) {
            chomp;
            next if(/(^\s*$)|(^\#)/); # skip blank or comment lines
            $count ++;
        }
        close $fh_in;
        return $count;
    }else {
        return 0;
    }
}

sub check_tools {
    my $tool = shift(@_);
    my %func = %{$tool};
    # @_ input the name of tool
    my @missing;
    for my $t (sort keys %func) {
        if( tool_path($t) ) {
            $func{$t} = tool_path($t);
        }else {
            push @missing, $t;
        }
    }
    if(@missing) {
        printf("The following tools are missing \n");
        for my $p (sort @missing) {
            printf("%-15s : Not found in \$PATH, %-30s\n\n", $p, '~/work/bin/temp/');
        }
    }
    my $st = (@missing)?0:1;
    return ($st, \%func);
#
    sub tool_path {
        my $tool = shift(@_);
        my $perldir = $ENV{HOME}. '/work/bin/temp';
        if(-e catfile($perldir, $tool)) {
            return catfile($perldir, $tool);
        }elsif(which($tool)) {
            return $tool;
        }else {
            return 0;
        }
    }
}

sub usage {
    die(qq/
Usage: parseBAM.pl <command> [arguments]\n
Command: stat   count mapping reads for each BAM
         view   create "*.bedgraph" files for each BAM
         tags   find tags from each BAM files (sRNA candidates)
\n/);
}

# Structure of the output directory:
# Output
#   |-stat_bams
#   |   |-bam.stat (output - 1)
#   |
#   |-view_bedgraph
#   |   |-sample1
#   |   |   |- *.fwd.bedgraph, *.rev.bedgraph, *.bam (output - 2)
#   |   |
#   |   |-sample2...
#   |
#   |-find_tags
#   |   |-sample1
#   |   |  |- sample1.tag.txt, *.tag.pos_sRNA.txt ....
#   |   |
#   |   |-sample2...
#   |   
#   |-merge_tags
#   |   |-Lib01
#   |   |  |-RNAz_out
#   |   |  |  |-best_RNAz.bed, best_hits.txt, SeqFA/
#   |   |  |-tag.newID.bed, tag.newID.pos_sRNA.count.txt, ...
#   |
#   |-report
#   |   |-Lib01.report.txt, Lib02.report.txt, ...       
#
#
# change log
# 1. convert bedgraph / hitogram graph using 'bedtools genomecov'
# 2. for strand-specific RNA-Seq, split PE BAM into two files: fwd.bam and rev.bam
# 3. find cov-tags with para: -cut-off: 100,
# 4. Using 6 MTB Complex genomes for RNAz analysis:
#      NC_008769: Mycobacterium bovis BCG str. Pasteur 1173P2
#      NC_015848: Mycobacterium canettii CIPT 140010059
#      NC_008596: Mycobacterium smegmatis str. MC2 155
#      NC_009525: Mycobacterium tuberculosis H37Ra
#      NC_000962: Mycobacterium tuberculosis H37Rv
#      NC_012943: Mycobacterium tuberculosis KZN 1435
# 5. set RNAz --cut-off=0.5
# 6. split these program into 3, for different purpose
#      stat: calculate total mapped read
#      view: create bedgraph files for genome browsers (eg: artemis, IGB)
#      tags: find cov regions in BAM files.
# 7. RNAz output: only report one region with the highest z-score if the input seq 
#    contain multiple regions with z-score > 0.5.
# 8. find cov regions (sRNA candidates): candidate should be 60 bp and 100 bp to its 
#    neighbor CDSs (genes).
# 9. you can input a custom GFF files (eg: only contain CDSs), to find sRNA candidates.
# 10. wrap count, TPM and RNAz score to one file in report directory.
# 
# 2014-11-03
# v0.1
#     1. support multiple bam files 
#     2. ONLY support single-chromosome sample
#     3. merge_tags by bedtools mrege: -s -d -1 -c 4,5,6 -o distinct,distinct,distinct
#     4. output results in out.dir/03.seqs
#
# 2015-02-07
# v0.2
#     1. support single-bam input file, Named: LibN
#     2. delete the step: copy reference file to current dir. insteat read the original fa/gff files
#
# 2015-04-05
# v0.3
#     1. Using 'HTSeq-count' instead of "bedtools multicov" to count reads on each features
#     2. Splite the BAM file into strand-specific files: fwd.bam and rev.bam
#        SE: 
#            samtools view -b -F 16 in.bam -o fwd.bam
#            samtools view -b -f 16 in.bam -o rev.bam
#        PE:
#            samtools view -b -f 128 -F 16 in.bam -o fwd1.bam
#            samtools view -b -f 80 in.bam -o fwd2.bam
#            samtools index fwd1.bam
#            samtools index fwd2.bam
#            samtools merge fwd.bam fwd1.bam fwd2.bam
#            #
#            samtools view -b -f 144 in.bam -o rev1.bam
#            samtools view -b -f 64 -F 16 in.bam -o rev2.bam
#            samtools index rev1.bam
#            samtools index rev2.bam
#            samtools merge rev.bam rev1.bam rev2.bam
#    3. Perform multiple sequence alignment by: ClustalW2 version 2.1, with default parameter
#    4. Perform RNAz analysis (RNAz version 2.1) --cut-off=0.5
#
# 2015-06-02
# v0.4
#    1. count reads on features by HTSeq-count: (find online note:)
#        SE: --stranded=yes
#        PE: --stranded=reverse (dUTP ssRNA-Seq)
#
# bug:
#    Line188: if input txt is blank, skip the 'paste' step, and copy the _sRNA to _count file. [BUG, not fix ye]
#
# 2015-06-15
# v0.5
#    1. replace HTSeq-count by FeatureCounts to count reads on each sRNAs. (much faster)
#
# Author: Wang Ming, wangmcas@gmail.com
