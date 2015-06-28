### !---------------------------------
### replaced by: chk_parseBAM.pl

#!/usr/bin/perl -w
use strict;
use warnings;

#######################################
# Extract cov regions from BAM file
# needs samtools and bedtools in $PATH
#######################################

use File::Basename qw(basename dirname);
use File::Spec::Functions qw(catfile catdir);
use File::Path qw(make_path remove_tree);
use File::Which;
use Cwd qw(abs_path cwd);
use POSIX qw(strftime);
use Getopt::Std;
use Data::Dumper;

### Global paras
my $cov_cutoff = 100;
my $ref;
my $gff;
my $outdir;
my $out_refseq;
#my $ref;
#my $gff;
my $out_tag;
my $out_merge;
my $out_seq;
my @bams;
my $bams_line;
my @bams_mapping;
### Tools
my $sort2bed;
my $cov2tag;
my $blast_to_rnaz;
my $sort2pos;
my $sort2cand;
my $tb_db = '/home/wangming/work/database/H37Rv/SixRv.fa';

&bam_to_tags($ref, $gff, @bams);

sub parse_para{
    my %opts = ();
    getopts("g:f:o:d:", \%opts);
    my $usage = "Usage: $0 [-f] <ref.fa> [-g] <ref.gff> [-o] <outdir> bam.list";
    die("$usage\n") if(@ARGV < 1);
    die "[-f] reference fasta not found.\n" unless defined  $opts{f};
    die "[-g] annotation GFF not found.\n" unless defined $opts{g};
    $tb_db = $opts{d} if defined $opts{d};
    $opts{o} = "search_tags" unless (defined $opts{o});
    $ref = $opts{f};
    $gff = $opts{g};
    $outdir = $opts{o};
    my $bam_list = shift(@ARGV);
    $bams_line = '';
    open my $list, "< $bam_list" or die "$!";
    while(<$list>){
        chomp;
        die "[$_] bam file not exists." unless -e $_;
        push @bams, abs_path($_);
        $bams_line .= ' '. abs_path($_);
    }
    close $list;
    @bams = sort @bams;
    check_tools();
    prep_wkdir();
}

sub check_tools{
    ### samtools and bedtools in $PATH
    die "[samtools] not found in the \$PATH" unless (which('samtools'));
    die "[bedtools] not found in the \$PATH" unless (which('bedtools'));
    ### Perl scripts
    $sort2bed = '/home/wangming/work/bin/sort2bed.pl';
    $cov2tag  = '/home/wangming/work/bin/temp/search_cov_regions.pl';
    $blast_to_rnaz = '/home/wangming/work/bin/temp/seq2rnaz.pl';
    $sort2pos = '/home/wangming/work/bin/temp/sort_to_position_sig.pl';
    $sort2cand= '/home/wangming/work/bin/temp/sort2candi_v1.pl';
#    $tb_db    = '/home/wangming/work/database/H37Rv/SixRv.fa';
    foreach my $p ($sort2bed, $cov2tag, $blast_to_rnaz, $sort2pos, $sort2cand, $tb_db){
        die "[$p] perl script not found." unless -e $p;
    }
}

sub prep_wkdir{
#    $ref = abs_path($ref);
#    $gff = abs_path($gff);
    system "samtools faidx $ref";
    # preparing wk dir
    $out_tag   = catdir($outdir, "01.tags");
    $out_merge = catdir($outdir, "02.merged");
    $out_seq   = catdir($outdir, "03.seqs");
    make_path($out_tag) unless -d $out_tag;
    make_path($out_merge) unless -d $out_merge;
    make_path($out_seq) unless -d $out_seq;
}

sub bam_to_tags{
    parse_para();
    @bams_mapping = &stat_bam(@bams);
    # 1 bam to tags
    print '['. &show_date . '] bam to tags' . "\n";
    my @tag_files;
    foreach my $b (sort @bams){
        my $tag = &parsing_bam($ref, $gff, $b);
        push @tag_files, $tag;
    }
    # 2. merge tags
    print '['. &show_date . '] merge tags' . "\n";
#    unlink glob "$out_merge/*.bed";
    my $merge_bed = &merge_tags(@tag_files);
    # 3. split merged file by type
    print '['. &show_date . '] split files' . "\n";
    my @lib_tags = &filter_seqs($merge_bed);
    my @lib_seqs = ();
    foreach my $g (@lib_tags){
        my $p = &tag_to_position($g);
        push @lib_seqs, $p;
    }
    # 4. Count TPM for each txt
    print '['. &show_date . '] count TPM' . "\n";
    my @tag_counts = &count_seq(@lib_seqs);
    my @tag_tpms = ();
    foreach my $c (@tag_counts){
        my $lib_tpm = &count_to_tpm($c);
        push @tag_tpms, $lib_tpm;
    }
    # 5. RNAz evaluation
    print '['. &show_date . '] run RNAz' . "\n";
    foreach my $f (@tag_tpms){
        my $tag_dir = dirname($f);
        my $rnaz_dir = catdir($tag_dir, "RNAz_out");
        make_path($rnaz_dir) unless -d $rnaz_dir;
        my $rnaz_file = &seq_to_rnaz($f, $rnaz_dir);
        my $comb_file = catfile($out_seq, basename($tag_dir). ".tpm_rnaz.txt" );
        my $tmp = &combine_rnaz_tpm($f, $rnaz_file, $comb_file);
    }
    # 6. Clear work path and Report
    
    print '['. &show_date . '] Finish'. "\n";
}

sub stat_bam{
    my @inputs = @_;
    my @stats = ();
    foreach my $b (sort @bams){
        my $mapped = qx(samtools idxstats $b | head -n1 |awk '{print \$3}');
        push @stats, $mapped;
    }
    return @stats;
}

sub parsing_bam{
    my ($ref, $gff, $bam) = @_;
    my ($prefix) = basename($bam) =~ /(.*)\.s\.bam/;
    die "[$bam] name is not *.s.bam" unless (defined $prefix);
    my $cov_n = catfile($out_tag, "$prefix\_coverage.n");
    my $cov_p = catfile($out_tag, "$prefix\_coverage.p");
    my $tag_n = catfile($out_tag, "$prefix\_tag.n");
    my $tag_p = catfile($out_tag, "$prefix\_tag.p");
    my $tag   = catfile($out_tag, "$prefix\_tag.txt");
    my $ref_fai = $ref . ".fai";
    my @sys_runs = ();
    push @sys_runs, "bedtools genomecov -ibam $bam -g $ref_fai -d -strand - > $cov_n";
    push @sys_runs, "perl $cov2tag -c $cov_cutoff -s - $cov_n > $tag_n";
    push @sys_runs, "bedtools genomecov -ibam $bam -g $ref_fai -d -strand + > $cov_p";
    push @sys_runs, "perl $cov2tag -c $cov_cutoff -s + $cov_p > $tag_p";
    push @sys_runs, "cat $tag_n $tag_p > $tag";
    foreach my $run (@sys_runs){
        system "$run";
    }
    return $tag;
}

sub merge_tags{
    my @txts = @_;
    my @beds = ();
    my $count = 1;
    my @sys_runs = ();
    my $tmp = catfile($outdir, "tmp");
    foreach my $m (@txts){
        # add prefix
        if(@txts == 1){
            push @sys_runs, "sed -e 's/^/LibN\_/' $m > $tmp";
        }else{
            push @sys_runs, "sed -e 's/^/Lib0$count\_/' $m > $tmp";
        }
        $count ++;
        my $m_bed = basename($m);
        $m_bed =~ s/\.txt/\.bed/;
        $m_bed = catfile($out_merge, $m_bed);
        push @sys_runs, "perl $sort2bed -t sort2bed -i $tmp -o $m_bed";
        push @beds, $m_bed;
        unlink $tmp;
    }
    my $bed_line = join " ", @beds;
    my $bed_all = catfile($out_merge, "all.sorted.bed");
    my $bed_merge = catfile($out_merge, "merged.bed");
    push @sys_runs, "cat $bed_line | sort -k1,1 -k2,2n > $bed_all";
    push @sys_runs, "bedtools merge -s -d -1 -c 4,5,6 -o distinct,distinct,distinct -i $bed_all > $bed_merge";
    foreach my $run (@sys_runs){
        system "$run";
    }
    return $bed_merge;
}

sub filter_seqs{
    my $in = shift(@_);
    my %lib;
    my %seq;
    my $counter = 1;
    foreach (@bams){
        my $lib_dir;
        if(@bams == 1){
            $lib_dir = catdir($out_merge, 'LibN');
        }else{
            $lib_dir = catdir($out_merge, "Lib0$counter");
        }
        make_path($lib_dir) unless -d $lib_dir;
        $lib{'tag'}->{$counter} = catfile($lib_dir, "tag.bed");
        $lib{'del'}->{$counter} = catfile($lib_dir, "del.bed");
        $counter ++;
    }
    $lib{'other'}->{0} = catfile($out_merge, 'unfiltered.txt');
    @{$seq{'other'}->{0}} = ();
    open F, "< $in" or die "$!";
    while(<F>){
        my $line = $_;
        my ($start, $end) = (split /\t/, $_)[1,2];
        my $length = $end - $start + 1;
        if(/Lib04/){
            if($length >=140){
                push @{$seq{'tag'}->{4}}, $line;
            }else{
                push @{$seq{'del'}->{4}}, $line;
            } 
        }elsif(/Lib03/){
            if($length >= 40 && $length <= 140){ 
                push @{$seq{'tag'}->{3}}, $line;
            }else{
                push @{$seq{'del'}->{3}}, $line;
            }
        }elsif(/Lib02/){
            if($length >= 40 && $length <= 80){
                push @{$seq{'tag'}->{2}}, $line;
            }else{
                push @{$seq{'del'}->{2}}, $line;
            }
        }elsif(/Lib01/){
            if($length >= 20 && $length <= 40){
                push @{$seq{'tag'}->{1}}, $line;
            }else{
                push @{$seq{'del'}->{1}}, $line;
            }
        }elsif(/LibN/){
            if($length >= 20){
                push @{$seq{'tag'}->{1}}, $line;
            }else{
                push @{$seq{'del'}->{1}}, $line;
            }
        }else{
            push @{$seq{'other'}->{0}}, $line;
        }
    }
    close F;
    foreach my $t (keys %lib){
        foreach my $n (keys %{$lib{$t}}){
            open OUT, "> $lib{$t}->{$n}" or die "$!";
            if(exists $seq{$t}->{$n}){
                print OUT join "", @{$seq{$t}->{$n}};
            }
            close OUT;
        }
    }
    return (sort values %{$lib{'tag'}});
}

sub count_seq{
    my @inputs = @_;
    my @f_counts = ();
    my @sys_runs = ();
    foreach my $f (@inputs){
        my $f_count = $f;
        $f_count =~ s/\.txt/\_count\.txt/;
        my $b1 = catfile(dirname($f), "tmp1.bed");
        my $b2 = catfile(dirname($f), "tmp2.bed");
        push @sys_runs, "perl $sort2bed -t sort2bed -i $f -o $b1";
        push @sys_runs, "bedtools multicov -s -bams $bams_line -bed $b1 > $b2";
        push @sys_runs, "perl $sort2bed -t bed2sort -i $b2 -o $f_count";
        push @f_counts, $f_count;
    }
    foreach my $run (@sys_runs){
        system "$run";
    }
    return @f_counts;
}

sub count_to_tpm{
    my $in_count = shift(@_);
    my $out_tpm = $in_count;
    $out_tpm =~ s/\_count\.txt/\_tpm.txt/;    
    open my $fh_count, "< $in_count" or die "$!";
    open my $fh_tpm, "> $out_tpm" or die "$!";
    while(<$fh_count>){
        chomp;
        my @tabs = split "\t";
        my @tpms = ();
        for(my $i = 0; $i < @bams; $i++){
            $tpms[$i] = sprintf "%.2f", $tabs[$i+12] * 1e6 / $bams_mapping[$i];
        }
        my $tag_id = shift(@tabs);
        my $tag_new_id = $tag_id;
        $tag_new_id = (split /\,/, $tag_id)[0] if($tag_id =~ /\,/); # choose the first id
        print $fh_tpm join "\t", ($tag_new_id, @tabs, @tpms, $tag_id);
        print $fh_tpm "\n";
    }
    close $fh_count;
    close $fh_tpm;
    return $out_tpm;
}

sub tag_to_position{
    my $in_bed = shift(@_);
    my $f1 = my $f2 = my $f3 = $in_bed;
    $f1 =~ s/\.bed/\.txt/;
    $f2 =~ s/\.bed/\.pos.txt/;
    $f3 =~ s/\.bed/\.pos_sRNA\.txt/;
    my @sys_runs = ();
    push @sys_runs, "perl $sort2bed -t bed2sort -i $in_bed -o $f1";
    push @sys_runs, "perl $sort2pos -f $gff -g $ref $f1 > $f2";
    push @sys_runs, "perl $sort2cand $f2";
    foreach my $run (@sys_runs){
        system($run);
    }
    return $f3;
}

sub seq_to_rnaz{
    my ($in, $rnaz_dir) = @_;
    my $in_fa = $in;
    $in_fa =~ s/\.txt/\.fa/;
    my $rnaz_output = $rnaz_dir;
    my $rnaz_log = catfile($rnaz_dir, "rnaz.log");
    my $rnaz_file = catfile($rnaz_dir, "best_RNAz.bed");
    my @sys_runs = ();
    push @sys_runs, "perl $sort2bed -t sort2fa -i $in -o $in_fa -g $ref";
    push @sys_runs, "perl $blast_to_rnaz -d $tb_db -o $rnaz_dir $in_fa > $rnaz_log ". '2>&1';
    foreach my $run (@sys_runs){
        system "$run";
    }
    return $rnaz_file;
}

sub combine_rnaz_tpm{
    my ($tpm, $rnaz_file, $seq_out) = @_;
    my %tpm;
    my %rnaz;
    %tpm = &read_txt($tpm);
    %rnaz = &read_rnaz($rnaz_file);
    open my $fh_cm, "> $seq_out" or die "$!";
    my $z_score = '-';
    foreach my $id (keys %tpm){
         my @tabs = split /\t/, $tpm{$id};
         my $old_id = pop @tabs;
         my $z_score = '-';
         $z_score = $rnaz{$id} if(exists $rnaz{$id});
         my $new_line = join "\t", ($id, @tabs, $z_score, $old_id);
         print $fh_cm $new_line, "\n";
    }
    close $fh_cm;
    return $seq_out;
}

sub read_txt{
    my $in = shift(@_);
    my %sort = ();
    open F, "<$in" or die "$!";
    while(<F>){
        chomp;
        next if(/^\s*$/); # skip blank lines
        my @tabs = (split /\t/, $_);
        my $id = shift(@tabs);
        $sort{$id} = join "\t", (@tabs);
    }
    close F;
    return %sort;
}

sub read_rnaz{
    my $in = shift(@_);
    my %rnaz = ();
    open my $fh, "<$in" or die "$!";
    while(<$fh>){
        chomp;
        next if(/^\s*$/); # skip blank lines
        my ($id, $score) = (split /\t/, $_)[0,4];
        my $id_new = $id;
        $id_new =~ s/^\>[a-zA-Z0-9]+\_//;
        $rnaz{$id_new} = $score;
    }
    close $fh;
    return %rnaz;
}

sub show_date{
    my $date = strftime "%Y-%m-%d %H:%M:%S", localtime;
    return $date;
}

__END__

Usage: bam_to_tags.pl -f ref.fa -g ref.gff -o out.dir bam.list

Output:

  out.dir
    |--01.tags
    |    |--*.coverage.n, *.coverage.p
    |    |--*.tag.n, *.tag.p, *.tag.txt
    |
    |--02.merged
    |    |--*merged.bed, all.sorted.txt, unfiltered.txt
    |    |--LibN[01-0N]
    |        |--tag*.sRNA.txt, tag*count.txt, tag*tpm.txt, ...
    |        |--RNAz_out
    |             |--best_RNAz.out
    |             |--best_hits.txt
    |
    |--03.seqs
    |    |--LibN.tpm.rnaz.txt
    |
    |--refseq
    |    |--ref.fa, ref.fa.fai, ref.gff 

Wang Ming wangmcas@gmail.com

Version:

v0.1
    1. support find cov regions from bam files [2-N].
    2. support for 1-chromosome strain.
    3. merge_tags by bedtools mrege: -s -d -1 -c 4,5,6 -o distinct,distinct,distinct 
    4. output results in out.dir/03.seqs
    
v0.2
    1. support 1-bam input, assign lib name: LibN
    2. delete the step: copy "refseq", use the original fa/gff.

v0.3
    1. replace bedtools count by "HTSeq-count"

