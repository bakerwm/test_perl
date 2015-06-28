### !-----------------------------
### replaced by: chk_parseBAM.pl

#!/usr/bin/perl -w
use strict;
use warnings;

use File::Basename qw(basename dirname);
use File::Spec::Functions qw(catfile catdir);
use File::Path qw(make_path remove_tree);
use Cwd qw(abs_path cwd);
use Getopt::Std;
use POSIX qw(strftime);

my %opts = ();
getopts("f:e:b:m:", \%opts);
die "Usage: perl $0 [-f] <genome.dir> [-e] <exclude.gene> [-b] <bam.dir> [-m] <Mapping.stat> <cov_dir> file1 file2 ... " if(@ARGV < 1);

my $cutoff_cov = 100;
my $ref_dir = abs_path($opts{f});
my $bam_dir = abs_path($opts{b});
my $map_stat = abs_path($opts{m});
my $exclude_genes = abs_path($opts{e});
my $ref_gff = '';
my $ref_fa = '';

my $cov_dir = shift;
my @cov_files = &parse_para;
my %f_perl = &find_perl();
my $cmd_line = join " ", ("perl", $0, '-f', $opts{f}, '-e', $opts{e}, '-b', $opts{b}, '-m', $opts{m}, $cov_dir);

### main script ###
my $outdir = "tags_out";
mkdir $outdir unless -d $outdir;
my $run_log = catfile($outdir, "run.log");

open my $fh, "> $run_log" or die "$!";
print $fh "Search cov_regions and count reads\nCommand line is:\n";
print $fh $cmd_line, "\n\n";

# find_sRNAs
my @tag_files = ();
foreach my $cov_n (sort @cov_files) {
    my $cov_p = $cov_n;
    $cov_p =~ s/\.n/.p/;
    my $fname = basename($cov_n);
    $fname =~ s/\.coverage\.n//;
    my $cov_tag_n = catfile($outdir, "$fname\_tag.n");
    my $cov_tag_p = catfile($outdir, "$fname\_tag.p");
    my $cov_tag_txt = catfile($outdir, "$fname\_tag.txt");
    print $fh '['. &show_date() . '] Search tags for '. $fname. "\n";
    &cov_to_regions($cov_n, $cov_tag_n);
    &cov_to_regions($cov_p, $cov_tag_p);
    system "cat $cov_tag_n $cov_tag_p > $cov_tag_txt";
    push @tag_files, basename($cov_tag_txt);
}
print $fh '['. &show_date() . '] Merge tags and count TPM '. "\n";

# tags to count
my @mapped = &check_mapping_stat($opts{m});

chdir $outdir or die "Cannot chdir to $outdir: $!";
my $tag_log = "tags_to_count.log";
&tags_to_count($tag_log, @tag_files);

print $fh '['. &show_date() . '] Run RNAz analysis.' . "\n";
system "sh run_rnaz.sh  > run_rnaz.sh.log";
system "find ./0* -name \"*TPM.bed\" > TPM.list ";
&combine_TPM("TPM.list");

# run RANz
chdir "../" or die "$!";
print $fh '['. &show_date() . '] Finish!'. "\n";
close $fh;

### Subroutines ###
sub parse_para {
    my $ref_dir = $opts{f};
    die "[$ref_dir] ref_dir not found." unless -d $ref_dir;
    my @gffs = glob"$ref_dir/*.gff";
    die "Make sure one [*gff] in [$ref_dir]" if(@gffs != 1);
    $ref_gff = shift(@gffs);
    my @fas = glob"$ref_dir/*.fna";
    die "Make sure one [*fna] in [$ref_dir]" if(@fas != 1);
    $ref_fa = shift(@fas);
    die "[$opts{e}] exclude_gene file not found." unless -e $opts{e};
    die "[$opts{b}] bam_dir not found." unless -d $opts{b};
    my @bams = <$opts{b}/*.s.bam>;
    die "[*.s.bam] bam_file not found." if (@bams < 1);
    die "[$opts{m}] mappint_stat file not found." unless -e $opts{m};
    my @cov_files = glob"$cov_dir/*coverage.n";
    die "[$cov_dir] no *.coverage.n files found" if(@cov_files < 1);
    return @cov_files;
}

sub show_date {
    my $date = strftime "%Y-%m-%d %H:%M:%S", localtime;
    return $date;
}

sub cov_to_regions {
    my ($in, $out) = @_;
    my $cov_cutoff = 100;
    my $strand = '';
    $strand = '+' if($in =~ /\.p$/);
    $strand = '-' if($in =~ /\.n$/);
    system "perl $f_perl{'search_cov_regions'} -c $cov_cutoff -s $strand $in > $out";
}

sub txt_to_position {
    my ($in, $out) = @_;
    system "perl $f_perl{'sort_to_position'} -f $ref_gff -e $exclude_genes $in > $out";
}

sub count_tags {
    my @ins = @_;
    my $in_list = join " ", @ins;    
    system "perl $f_perl{'tags_to_count'} -f $ref_dir -b $bam_dir -m $map_stat -c no $in_list";
}

sub combine_TPM {
    my $in = shift(@_);
    system "perl $f_perl{'combine_TPM'} $in"
}

sub find_perl {
    my $p1 = '/home/wangming/work/bin/temp/search_cov_regions.pl';
    my $p2 = '/home/wangming/work/bin/get_sRNA/sort_to_position_sig.pl';
    my $p3 = '/home/wangming/work/bin/get_sRNA/sort2candi_v1.pl';
    my $p4 = '/home/wangming/work/bin/sort2bed.pl';
    my $p5 = '/home/wangming/work/bin/temp/combine_TPMs.pl';
    my $p6 = '/home/wangming/work/bin/temp/count_to_TPM.pl';
    my $p7 = '/home/wangming/work/bin/temp/blast2RNAz.pl';
    my $p8 = '/home/wangming/work/database/H37Rv/SixRv.fa';

    my %find_perl = ('search_cov_regions' => $p1,
                     'sort_to_position' => $p2,
                     'sort_to_candi' => $p3,
                     'sort_to_bed' => $p4,
                     'combine_TPM' => $p5,
                     'count_to_TPM' => $p6,
                     'blast_to_RNAz' => $p7,
                     'SixRv' => $p8);
    foreach my $p (keys %find_perl) {
        die "[$find_perl{$p}] not found." unless -e $find_perl{$p};
    }
    return %find_perl;
}

### BEGIN tags_to_count ###
sub tags_to_count {
    my $tags_to_count_log = shift(@_);
    my @infiles = @_;
    open my $tag_log, "> $tags_to_count_log" or die "$!";
    ## clear wrok dir ##
    my $clear_wk = '# Clear workspace: 01.merge_no/RNAz_out/* and 02.merge_yes/RNAz_out/*' . 
                   "\n" . 'rm -r 0*/RNAz_out/*' . "\n\n";
    while(my $n = <0*merge*/Lib*>) {
        remove_tree("$n");
    }
    ## working ##
    print $tag_log '['. &show_date() . '] not_merged tags: count'. "\n";
    my $merge_type = '01.merge_no';
    my $cmd1 = &files_to_TPM($merge_type, @infiles);
    print $tag_log '['. &show_date() . '] merged tags: count', "\n";
    $merge_type = '02.merge_yes';
    my @wk_files = &merge_files_bed($merge_type, @infiles);
    my $cmd2 = &files_to_TPM($merge_type, @wk_files);
    open OUT, "> run_rnaz.sh" or die "$!";
    print OUT join "\n",($clear_wk, $cmd1, $cmd2), "\n";
    close OUT;
    print $tag_log '[' . &show_date() . '] Finish!', "\n";
    close $tag_log;
}

sub check_mapping_stat{
    my $map_ex = "Name\tTotal\tMapped\tRatio\n18-40\t1000\t100\t0.1\n";
    my $m = shift(@_);
    my @maps = ();
    open my $map, "< $m" or die "$!";
    while(<$map>) {
        chomp;
        next unless(/^\d+/);
        my @tmp = split /\s+/;
        die "Example mapping file:\n$map_ex\n" unless(@tmp == 4 && $tmp[3] =~ /^\d+$/);
        push @maps, $tmp[2];
    }
    return @maps;
}

sub merge_files_bed{ # merge different txt files, then split by length
    my $merge_type = shift(@_);
    make_path($merge_type);
    my @infiles = @_;
    warn "No need to merge for one file" if(@infiles < 2);
    foreach my $f (@infiles) {
        my $t = basename($f);
        &add_id_prefix($f, "$t\.id");
    }
    system "cat *.id > All_Libs_id.txt";
    &txt_to_bed("All_Libs_id.txt");
    system "sort -k1,1 -k2,2n All_Libs_id.bed -o All_Libs_id.bed";
    my $merged_file = catfile('02.merge_yes', 'merged.bed');
    system "bedtools merge -s -d -1 -c 4,5,6 -o distinct,distinct,distinct -i All_Libs_id.bed > $merged_file";
    &bed_to_txt($merged_file);
    system "rm  *.id  All_Libs_id.*";
    my ($strain) = basename($infiles[0]) =~ /^(\w+)\./;
    my $merged_txt = $merged_file;
    $merged_txt =~ s/\.bed/.txt/;
    my @can_files = &split_merge_by_type($merged_txt, $merge_type, $strain);
    return @can_files;
}

sub files_to_TPM {
    my $merge_type = shift(@_);
    my @infiles = @_;
    my @cmd_lines = ();
    foreach my $i (@infiles) {
        my $chk = &check_file($i);
        next unless($chk); # have non-zero line
        my ($f) = basename($i) =~ /(^\w+\.\d+)/;
        my ($lib, $lib_range, $bam, $total_reads) = &parse_lib_bam($i);
        my $outpath = catdir($merge_type, $lib);
        make_path($outpath) unless -d $outpath;
        my $f_pre = catfile($outpath, $f);
        system "cp $i  $f_pre\-0.txt";
        &size_filter("$f_pre\-0.txt", "$f_pre\-1.txt", $lib_range);
        &txt_to_position("$f_pre\-1.txt", "$f_pre\-2.txt");
        &txt_to_candi("$f_pre\-2.txt");
        system "cp -r $f_pre\-2_sRNA.txt $f_pre\-3.txt";
        &txt_to_bed("$f_pre\-3.txt");
        system "cut -f1-6 $f_pre\-3.bed > t.bed && mv t.bed $f_pre\-3.bed";
        system "bedtools multicov -s -bams $bam -bed $f_pre\-3.bed > $f_pre\-3_count.bed";
        system "perl $f_perl{count_to_TPM} -i $f_pre\-3_count.bed -n $total_reads > $f_pre\-3_TPM.bed";
        system "perl $f_perl{sort_to_bed} -t bed2sort -i $f_pre\-3_count.bed -o $f_pre\-3_count.txt";
        system "perl $f_perl{sort_to_bed} -t sort2fa  -i $f_pre\-3_count.txt -o $f_pre\-3_count.fa -g $ref_fa";
        push @cmd_lines, "\n".'# Perform RNAz analysis';
        my $rnaz_path = catdir($outpath, "RNAz_out");
        make_path($rnaz_path) unless -d $rnaz_path;
        push @cmd_lines, "perl $f_perl{blast_to_RNAz} -d $f_perl{SixRv} -o $rnaz_path $f_pre\-3_count.fa > rnaz.log \n";
    }
    return join "\n", @cmd_lines;
}

sub parse_lib_bam{
    my $in = shift(@_);
    my $f_name = basename($in);
    my ($lib, $lib_range, $total_reads, $bam);
    my @tmps;
    if ($f_name =~ /\.01/) {
        $lib = "Lib01_18-40";
        $lib_range = "18:40";
        $total_reads = $mapped[0];
        @tmps = <$bam_dir\/*01.s.bam>;
        $bam = $tmps[0];
    }elsif ($f_name =~ /\.02/ ) {
        $lib = "Lib02_40-80";
        $lib_range = "40:80";
        $total_reads = $mapped[1];
        @tmps = <$bam_dir\/*02.s.bam>;
        $bam = $tmps[0];
    }elsif ($f_name =~ /\.03/) {
        $lib = "Lib03_80-140";
        $lib_range = "80:140";
        $total_reads = $mapped[2];
        @tmps = <$bam_dir\/*03.s.bam>;
        $bam = $tmps[0];
    }elsif ($f_name =~ /\.04/) {
        $lib = "Lib04_140";
        $lib_range = "140:1000000";
        $total_reads = $mapped[3];
        @tmps = <$bam_dir\/*04.s.bam>;
        $bam = $tmps[0];
    }else {
    }
    return ($lib, $lib_range, $bam, $total_reads);
}

sub split_merge_by_type{ # split a file by the lenght
    my ($in, $merge_type, $strain) = @_; 
    # Clear the work dir
    while (my $n = <$merge_type\/Lib*>) {
        remove_tree($n);
    }
    # Lib01
    make_path("$merge_type\/Lib01_18-40") unless -d "$merge_type\/Lib01_18-40";
    system "grep \'LIB01\' $in | awk \'{if(\$3 < 40) print}\' > $merge_type\/Lib01_18-40\/$strain.01-00.txt ";
    # Lib02
    make_path("$merge_type\/Lib02_40-80") unless -d "$merge_type\/Lib02_40-80";
    system "grep \'LIB02\' $in | awk \'{if(\$3 >= 40 && \$3 < 80) print}\' > $merge_type\/Lib02_40-80\/$strain.02-00.txt ";
    # Lib03
    make_path("$merge_type\/Lib03_80-140") unless -d "$merge_type\/Lib03_80-140"; 
    system "grep \'LIB03\' $in | awk \'{if(\$3 >= 80 && \$3 < 140) print}\' > $merge_type\/Lib03_80-140\/$strain.03-00.txt ";
    # Lib04
    make_path("$merge_type\/Lib04_140") unless -d "$merge_type\/Lib04_140";
    system "grep \'LIB04\' $in | awk \'{if(\$3 > 140) print}\' > $merge_type\/Lib04_140\/$strain.04-00.txt ";
    my @c_files = glob"$merge_type\/Lib*\/*.txt";
    return @c_files;
}

sub add_id_prefix{
    my ($in, $out) = @_;
    my $in_name = basename($in);
    if ($in_name =~ /\.01\_/){
        system "sed -e \'s/^/LIB01_/\' $in |cut -f1-6 > $out";
    } elsif ($in_name =~ /\.02\_/) {
        system "sed -e \'s/^/LIB02_/\' $in |cut -f1-6 > $out";
    } elsif ($in_name =~ /\.03\_/) {
        system "sed -e \'s/^/LIB03_/\' $in |cut -f1-6 > $out";
    } elsif ($in_name =~ /\.04\_/) {
        system "sed -e \'s/^/LIB04_/\' $in |cut -f1-6 > $out";
    } else {
    
    }
}

sub txt_to_bed{
    my $in = shift(@_);
    &check_file($in);
    my $out = $in;
    $out =~ s/\.txt/.bed/;
    system "perl $f_perl{sort_to_bed} -t sort2bed -i $in -o $out";
}

sub bed_to_txt{
    my $in = shift(@_);
    &check_file($in);
    my $out = $in;
    $out =~ s/\.bed/.txt/;
    system "perl $f_perl{sort_to_bed} -t bed2sort -i $in -o $out";
}

sub txt_to_candi{
    my ($in) = @_;
    &check_file($in);
    system "perl $f_perl{sort_to_candi} $in";
}

sub size_filter{ # select seqs by range of length
    my ($in, $out, $range) = @_; # in.txt out.txt 18-40: >140 "140:N", <40 "N:40"
    &check_file($in);
    my ($min, $max) = split "\:",$range;
    mkdir "tmp" unless -d "tmp";
    if($min =~ /^\d+$/){
        system "awk \'{if(\$3 >= $min) print}\' $in > tmp/temp01";
        if($max =~ /^\d+$/){
            system "awk \'{if(\$3 < $max) print}\' tmp/temp01 > $out";
        }else{
            system "mv tmp/temp01.txt $out";
        }
    }else{
        die "Need at least one NUM $range" unless($max =~ /^\d+$/);
        system "awk \'{if(\$3 < $min) print}\' $in > $out.txt";
    }
    remove_tree("tmp");
}

sub check_file{
    my $in = shift(@_);
    my $Length = 0;
    warn "[$in] not found" unless -e $in;
    my $line = `wc -l $in`;
    ($Length) = $line =~ /^(\d+)/;
    return($Length);
}
### END tags_to_count ###
