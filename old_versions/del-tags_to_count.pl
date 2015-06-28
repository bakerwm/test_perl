### !------------------------------------------------
### replaced by: chk_feature2count.pl

#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

use File::Basename qw(dirname basename);
use File::Spec::Functions qw(catdir catfile);
use File::Path qw(make_path remove_tree);

my %opts = ();
$opts{c} = 'no';
getopts("f:b:m:c:", \%opts);
die ("Usage: tags_to_count.pl [-f] <ref_dir> [-b] <bam_dir> [-m] <mapping.stat> [-c] <countTPM: yes|no> <file1> <file2> ...\n") if (@ARGV == 0);

#============ Start Pre-wk =============#
## Scripts
my $sort_to_candi = '/home/wangming/work/bin/get_sRNA/sort2candi_v1.pl';
my $sort_to_position = '/home/wangming/work/bin/get_sRNA/sort_to_position_sig.pl';
my $sort_to_bed = '/home/wangming/work/bin/sort2bed.pl';
my $count_to_TPM = '/home/wangming/work/bin/temp/count_to_TPM.pl';
my $blast_to_RNAz = '/home/wangming/work/bin/temp/blast2RNAz.pl';
my $combine_TPMs = '/home/wangming/work/bin/temp/combine_TPMs.pl';
my $rnaz_db = '/home/wangming/work/database/H37Rv/SixRv.fa';
## check script existence
my @scripts = ($sort_to_candi, $sort_to_position, $sort_to_bed, $count_to_TPM, $blast_to_RNAz, $combine_TPMs);
foreach my $s(@scripts){
   die "$s, $!" unless -e $s;
}

## Parse ref & BAMs
my $ref_dir = $opts{f};
my @ref_gffs = <$ref_dir\/*.gff>;
my @ref_fas = <$ref_dir\/*.fna>;
die "Need only one [*gff|*fna] file in [$ref_dir]" if(@ref_gffs != 1 || @ref_fas != 1);
my $ref_gff = shift(@ref_gffs);
my $ref_fa = shift(@ref_fas);
my $bam_dir = $opts{b};
my @bam_files = glob"$bam_dir\/*.s.bam";
die ("Cannot find BAM files [*.s.bam]: $bam_dir") if (scalar(@bam_files < 1));

## Mapping stat file
my @mapped = &check_mapping_stat($opts{m});

#============ END Pre-wk ================#

#============ START main script ===========#
## Prepare dir ##
my @cmds;
push @cmds, '# merge tags and count TPM';
push @cmds, '# The command line is: ';
push @cmds, join ' ', ('#', $0, '-f', $opts{f}, '-b', $opts{b}, '-m', $opts{m}, @ARGV);
push @cmds, '# Clear the old files and folders: *merge*/Lib*';
print "Clear wk dir\n";
while(my $n = <./*merge*/Lib*>) {
    remove_tree("$n");
}
push @cmds, "\n".'# For not merged files';
push @cmds, '# Transform input file from txt to bam format';
print 'Transform input file from txt to bam format', "\n";

## working ##
print 'For 01.merge_no folder' . "\n";
my $merge_type = '01.merge_no';
my $cmd1 = &files_to_TPM($merge_type, @ARGV);

print 'For 02.merge_yes folder' . "\n";
$merge_type = '02.merge_yes';
my $run_rnaz_sh = 'run_rnaz.sh';
my @wk_files = &merge_files_bed($merge_type, @ARGV);
my $cmd2 = &files_to_TPM($merge_type, @wk_files);
open OUT, "> $run_rnaz_sh" or die "$!";
push @cmds, ($cmd1,"\n", $cmd2, "\n");
print OUT join "\n", @cmds;
close OUT;

## Count reads TPM and RNAz analysis
if($opts{c} =~ /^yes|y$/i) {
   print 'Count reads, TPM and RNAz analysis' . "\n";
   system "sh $run_rnaz_sh  > $run_rnaz_sh\.log";
   system "find ./0* -name \"*TPM.bed\" > TPM.list ";
   system "perl $combine_TPMs TPM.list";
} else {
   print "Skip RNAz analysis\n[Run: sh $run_rnaz_sh > $run_rnaz_sh\.log]" . "\n";
}

#============ END main script ================#

#============ Subroutines ================#
sub check_mapping_stat{
    my $map_ex = "ID\tName\tTotal\tMapped\n1\t18-40\t1000\t100\n";
    my $m = shift(@_);
    my @maps = ();
    open my $map, "< $m" or die "$!";
    while(<$map>) {
        chomp;
        next unless(/^\d+/);
        my @tmp = split /\s+/;
        die "Example mapping file:\n$map_ex\n" unless(@tmp == 4 && $tmp[3] =~ /^\d+$/);
        push @maps, $tmp[3];
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
    my ($strain) = basename($ARGV[0]) =~ /^(\w+)\./;
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
        &txt_to_position("$f_pre\-1.txt", "$f_pre\-2.txt", $ref_gff);
        &txt_to_candi("$f_pre\-2.txt");
        system "cp -r $f_pre\-2_sRNA.txt $f_pre\-3.txt";
        &txt_to_bed("$f_pre\-3.txt");
        system "cut -f1-6 $f_pre\-3.bed > t.bed && mv t.bed $f_pre\-3.bed";
        system "bedtools multicov -s -bams $bam -bed $f_pre\-3.bed > $f_pre\-3_count.bed";
        system "perl $count_to_TPM -i $f_pre\-3_count.bed -n $total_reads > $f_pre\-3_TPM.bed";
        system "perl $sort_to_bed -t bed2sort -i $f_pre\-3_count.bed -o $f_pre\-3_count.txt";
        system "perl $sort_to_bed -t sort2fa  -i $f_pre\-3_count.txt -o $f_pre\-3_count.fa -g $ref_fa";
        push @cmd_lines, "\n".'# Perform RNAz analysis';
        my $rnaz_path = catdir($outpath, "RNAz_out");
        make_path($rnaz_path) unless -d $rnaz_path;
        push @cmd_lines, "perl $blast_to_RNAz -d $rnaz_db -o $rnaz_path $f_pre\-3_count.fa > $rnaz_path\/rnaz.log \n";
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
    system "perl $sort_to_bed -t sort2bed -i $in -o $out";
}

sub bed_to_txt{
    my $in = shift(@_);
    &check_file($in);
    my $out = $in;
    $out =~ s/\.bed/.txt/;
    system "perl $sort_to_bed -t bed2sort -i $in -o $out";
}

sub txt_to_position{
    my ($in, $out, $gff) = @_;
    &check_file($in);
    system "perl $sort_to_position -i $in -o $out -f $gff";
}

sub txt_to_candi{
    my ($in) = @_;
    &check_file($in);
    system "perl $sort_to_candi $in";
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
#============ END Subroutines ==============#
