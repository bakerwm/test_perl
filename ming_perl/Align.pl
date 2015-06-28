#!/usr/bin/perl -w
use strict;

use lib qw(/home/wangming/localperl/share/perl5 ); # local module lib

use warnings;
use Cwd qw(abs_path cwd);
use File::Basename qw(basename dirname);
use File::Which qw(which);
use File::Spec::Functions qw(catdir catfile);
use Getopt::Long;
use Pod::Usage;
# History:2014-11-23  updated: line-199, 245, 281; elsif($lib eq "03" || $lib eq "04")  => elsif($lib > "02")

my $Aligner=1;
my $reads='';
my $reference='';
my $name='';
my $suffix='';
my $help='';
my $man=0;

GetOptions('Aligner|p=i' => \$Aligner,
           'Reads|r=s' => \$reads,
           'reference|d=s' => \$reference,
           'name|n=s' => \$name,
           'suffix|x=s' => \$suffix,
           'help|h' => \$help,
           'man' => \$man
        ) or pod2usage(-verbose => 0);

pod2usage(1) if ($help);
pod2usage(-verbose => 2) if ($man);

## Prepare output dir
my $original_dir = abs_path(cwd);
my $output_suffix = ($suffix eq '')?"":"\_$suffix";
my $Alignment_dir = catdir($original_dir, "02.Alignment$output_suffix");
my $Coverage_dir  = catdir($original_dir, "03.Coverage$output_suffix");
my $Merge_dir     = catdir($original_dir, "04.Merge_tags$output_suffix");
my $Count_dir     = catdir($original_dir, "05.Count_reads$output_suffix");
#&check_dir($Alignment_dir, $Coverage_dir, $Merge_dir, $Count_dir);

## Check tools
my $soap        = which('soap');
my $bowtie      = which('bowtie');
my $bowtie2     = which('bowtie2');
my $samtools    = which('samtools');
my $bedtools    = which('bedtools');
my $htseq_count = which('htseq-count');
&Check_file_existence($soap, 
#                      $bowtie, 
                      $bowtie2, 
                      $samtools, 
                      $bedtools, 
                      $htseq_count);

my $my_bin         = catdir('/home/wangming', 'work/bin/get_sRNA');
my $get_sRNA_pl    = catfile($my_bin, 'get_sRNA.pl');
my $getPosition_pl = catfile($my_bin, 'getPosition_v4.pl');
my $merge4all_pl   = catfile($my_bin, 'merge4all.pl');
my $sort2candi_pl  = catfile($my_bin, 'sort2candi_v1.pl');
my $sort2temp_pl   = catfile($my_bin, 'sort2temp.pl');
my $sort2bed_pl    = catfile($my_bin, 'sort2bed.pl');
my $soapPE_merge_pl= catfile($my_bin, 'soapPE_merge.pl');
my $soap2sam_pl    = catfile($my_bin, 'soap2sam.pl');
my $statMapping_pl = catfile($my_bin, 'statMapping.pl');
&Check_file_existence($get_sRNA_pl, $getPosition_pl, $merge4all_pl, $sort2candi_pl,
                      $sort2temp_pl, $sort2bed_pl, $soapPE_merge_pl, $soap2sam_pl,
                      $statMapping_pl);

## Check input Aligner
pod2usage(-verbose => 1, -message => "$0: -p | --Aligner, $Aligner not in [1-3]. \n") 
          unless ($Aligner == 1 || $Aligner == 2 || $Aligner == 3);

## Check reads
my $reads_dir = $reads;
my @list_reads = <$reads_dir/*.f*>;
die "$0: -r, [$reads_dir], No fa/fastq files found. \n" if(@list_reads < 1);
my %check_reads = &check_reads(@list_reads);
if($check_reads{'flag'} eq 'NA'){
    print "Found reads file:\n";
    print join"\n", @list_reads, "\n";
    die 'The name of reads file not match [\w+\d\d].fastq | [\w+\d\d\_\d].fastq. ', "\n";
}

my %reads_info = %{$check_reads{'libs'}};
my @reads_libs = keys %reads_info;

sub check_reads{
    my $flag = '';
    my %lib = ();
    my %out = ();
    foreach my $f (@_){
        my $f_name = basename($f);
        $flag = 'NA' unless ($f_name =~ /^\w+\.fa$|\w+\.fastq$/);
        if ($f_name =~ /^\w+\d\d\.fa|^\w+\d\d\.fastq$/){
            my ($num) = $f_name =~ /^\w+(\d\d)\.fa/;
            push @{$lib{$num}}, $f_name;
        }elsif($f_name =~ /^\w+\d\d\_\d\.fa|^\w+\d\d\_\d\.fastq$/){
            my ($num) = $f_name =~ /^\w+(\d\d)\_\d\.fa/;
            push @{$lib{$num}}, $f_name;;
        }
    }
    $out{'flag'} = $flag;
    %{$out{'libs'}} = %lib;
    return %out;
}

## Check reference fa, gff, fai
$reference = abs_path($reference);
my $gff = $reference;
$gff =~ s/\.f.*$/.gff/; 
$gff = abs_path($gff);
die "$0: -d | --reference, [$reference] Not exists.\n" unless (-e $reference);
die "$0: -d | --reference, [$gff] Not found. \n" unless(-e $gff);
`$samtools faidx $reference` unless -e "$reference\.fai";

## Sample name
die "$0: -n | --name, [$name] should be [a-zA-Z0-9_]. \n" unless ($name =~ /^\w+$/);

## Start generate Shell commands
my @cmd_lines = ();

my $pTime = 'Date=$(date)
echo "$Date  ';


check_dir($Alignment_dir, $Coverage_dir, $Merge_dir, $Count_dir);

# Step 1. Filter reads
push @cmd_lines, '#!/bin/bash';
push @cmd_lines, "\n$pTime \t \# Step 1. Filter reads \"  \n";

# Step 2. Create Alignment command lines
push @cmd_lines, "$pTime \t \# Step 2. Mapping reads and Find sRNAs. \"  \n";

$Alignment_dir = abs_path($Alignment_dir);
# Aligner_para
my $soap_SE_para    = '-v 3 -r 2 -p 4 -M 4';
my $soap_PE_para    = '-v 3 -r 2 -p 4 -m 40 -x 300 -s 40 -l 200';
my $bowtie_SE_para  = '-v 3 --best --strata -S -m 100 -X 300 --chunkmbs 256 -p 4';
my $bowtie_PE_para  = '-v 3 --best --strata -S -m 100 -I 60 -X 300 --chunkmbs 256 -p 4';
my $bowtie2_SE_para = '-p 4';
my $bowtie2_PE_para = '-I 50 -X 200 -p 4';

foreach my $lib (sort @reads_libs){
	my $sam_file   = catfile($Alignment_dir, "$name.$lib.sam");
    push @cmd_lines, "\n$pTime \t \#\#\# Mapping reads: $name $lib \"  \n ";
# soap commands
	if($Aligner == 1){
        # Check index file
        my $soap_index = $reference.".index";
        my @check_index = <$soap_index.*>;
        my $soap_build = which('2bwt-builder'); # print join "\n", @check_index, "\n"; exit;
        system "$soap_build $reference" if (@check_index < 1);

        my $soap_out = catfile($Alignment_dir, "$name.$lib.soap");
        my $soap_PE = catfile($Alignment_dir, "$name.$lib.soapPE");
        my $soap_PESingle = catfile($Alignment_dir, "$name.$lib.soapSingle");
        my $soap_un = catfile($Alignment_dir, "$name.$lib.unsoap");
        # Soap for SE
        if($lib eq "01" || $lib eq "02"){
            my $SE_reads_name = $reads_info{$lib}[0];
            my $SE_reads = abs_path(catfile($reads_dir, $SE_reads_name));
            my $run_soap = join " ", ($soap, $soap_SE_para, '-D', $soap_index, 
                    '-a', $SE_reads, '-o', $soap_out, '-u', $soap_un);
            push @cmd_lines, $run_soap;
        # Soap for PE
        }elsif($lib >= "03"){
            my @PE_reads = @{$reads_info{$lib}};
            my $PE_reads_1 = abs_path(catfile($reads_dir, $PE_reads[0]));
            my $PE_reads_2 = abs_path(catfile($reads_dir, $PE_reads[1]));
            my $run_soap = join " ", ($soap, $soap_PE_para, '-D', $soap_index, 
                    '-a', $PE_reads_1, '-b', $PE_reads_2, '-o', $soap_out, 
                    '-2', $soap_PESingle, '-u', $soap_un);
            my $soap_PE2SE_fa  = catfile($Alignment_dir, "$name.$lib.fa");
            my $soap_PE2SE_del = catfile($Alignment_dir, "$name.$lib.del");
            my $run_soapPEmerge = join " ", ('perl', $soapPE_merge_pl, '-a', $reference, 
                    '-i', $soap_PE, '-o', $soap_out, '-c', $soap_PE2SE_fa, '-d',
                    $soap_PE2SE_del );
            push @cmd_lines, ($run_soap, $run_soapPEmerge);
        }
        my $run_soap2sam = join " ", ('perl', $soap2sam_pl, $soap_out, '>', $sam_file);
        push @cmd_lines, ($run_soap2sam);
# bowtie commands
	}elsif($Aligner == 2){ # bowtie
	    # Check index file
        my $reference_dir = dirname($reference);
        my $bowtie_index = catfile($reference_dir, $name);
        my @check_index  = <$bowtie_index\.*\.ebwt>;
        my $bowtie_build = which('bowtie-build'); &Check_file_existence($bowtie_build);
        system"$bowtie_build -q -f $reference $bowtie_index" if(@check_index < 1);

        my $run_bowtie = '';
        my $read_format = '-q';
        if($lib eq '01' || $lib eq '02'){
            my $SE_reads_name = $reads_info{$lib}[0];
            my $SE_reads = abs_path(catfile($reads_dir, $SE_reads_name));
            $read_format = '-f' if($SE_reads =~ /\.fa$/);
            $run_bowtie = join " ", ($bowtie, $bowtie_SE_para, $bowtie_index, $read_format,
                    $SE_reads, $sam_file);
        }elsif($lib >= '03'){
            my @PE_reads = @{$reads_info{$lib}};
            my $PE_reads_1 = abs_path(catfile($reads_dir, $PE_reads[0]));
            my $PE_reads_2 = abs_path(catfile($reads_dir, $PE_reads[1]));
            $read_format = '-f' if($PE_reads_1 =~ /\.fa$/);
            $run_bowtie = join " ", ($bowtie, $bowtie_PE_para, $bowtie_index, 
                    '-1', $PE_reads_1, '-2', $PE_reads_2, $sam_file);        
        }else{
            die "Check the type of reads: $lib.";
        }
        push @cmd_lines, ($run_bowtie);
# bowtie2 commands
	}elsif($Aligner == 3){ # bowtie2
        # Check index file
        my $reference_dir = dirname($reference);
        my $bowtie2_index = catfile($reference_dir, $name);
        my @check_index = <$bowtie2_index\.*\.bt2>;
        my $bowtie2_build = which('bowtie2-build'); &Check_file_existence($bowtie2_build);
        system"$bowtie2_build -q -f $reference $bowtie2_index" if(@check_index < 1);

        my $run_bowtie2 = '';
        my $read_format = '-q';
        if($lib eq '01' || $lib eq '02'){
            my $SE_reads_name = $reads_info{$lib}[0];
            my $SE_reads = abs_path(catfile($reads_dir, $SE_reads_name));
            $read_format = '-f' if($SE_reads =~ /\.fa$/);
            $run_bowtie2 = join " ", ($bowtie2, $bowtie2_SE_para, '-x', $bowtie2_index, 
                    $read_format, $SE_reads, '-S', $sam_file);
        }elsif($lib >= '03'){
            my @PE_reads = @{$reads_info{$lib}};
            my $PE_reads_1 = abs_path(catfile($reads_dir, $PE_reads[0]));
            my $PE_reads_2 = abs_path(catfile($reads_dir, $PE_reads[1]));
            $read_format = '-f' if($PE_reads_1 =~ /\.fa$/);
            $run_bowtie2 = join " ", ($bowtie2, $bowtie2_PE_para, '-x', $bowtie2_index,
                    '-1', $PE_reads_1, '-2', $PE_reads_2, '-S', $sam_file);        
        }else{
            die "Check the type of reads: $lib.";
        }
        push @cmd_lines, ($run_bowtie2);
# Other Aligner
	}else{
	    die "$0: -p | --Aligner [$Aligner] Not in [soap, bowtie, bowtie2] \n";	
	}
# Run sam to bam
    my $sam_prefix = catfile($Alignment_dir, "$name.$lib");
#    my $sam_file = $sam_prefix.'.sam';
    my $bam_file = $sam_prefix.'.bam';
    my $bam_s_file = $sam_prefix.'.s.bam';

    my $run_sam2bam  = join " ", ($samtools, 'view -bt', "$reference\.fai", $sam_file, 
    '-o', $bam_file);
    my $run_bam2sort = join " ", ($samtools, 'sort', $bam_file, "$sam_prefix\.s");
    my $run_bam2index = join " ", ($samtools, 'index', $bam_s_file);
    push @cmd_lines, ($run_sam2bam, $run_bam2sort, $run_bam2index);
# Run bam to cov, find RNA
    push @cmd_lines, "\n$pTime \t \#\#\# Find sRNA. $name $lib\"  \n ";
    my $Coverage_n = catfile($Coverage_dir, "$name.$lib.coverage.n");
    my $Coverage_p = catfile($Coverage_dir, "$name.$lib.coverage.p");
    my $tag_tmp_n  = catfile($Coverage_dir, "$name.$lib\_tag\.tmp\.n");
    my $tag_tmp_p  = catfile($Coverage_dir, "$name.$lib\_tag\.tmp\.p");
    my $tag_tmp    = catfile($Coverage_dir, "$name.$lib\_tag\.tmp");
    my $tag_txt    = catfile($Coverage_dir, "$name.$lib\_tag\.txt");
    my $run_bam2cov_n = join " ", ($bedtools, 'genomecov', '-ibam', $bam_s_file, '-d -g', 
    "$reference\.fai", '-strand - >', $Coverage_n);
    my $run_bam2cov_p = join " ", ($bedtools, 'genomecov', '-ibam', $bam_s_file, '-d -g', 
    "$reference\.fai", '-strand + >', $Coverage_p);
    my $run_cov2RNA_n = join " ", ('perl', $get_sRNA_pl, '-i', $Coverage_n, 
    '-c 0.5 -s - -o', $tag_tmp_n);
    my $run_cov2RNA_p = join " ", ('perl', $get_sRNA_pl, '-i', $Coverage_p,
    '-c 0.5 -s + -o', $tag_tmp_p);
    my $run_catNP = join " ", ('cat', $tag_tmp_n, $tag_tmp_p, '>', $tag_tmp);
    my $run_tag2pos = join " ", ('perl', $getPosition_pl, '-i', $tag_tmp, '-f', $gff, 
    '-o', $tag_txt);
    my $run_sort2candi = join " ", ('perl', $sort2candi_pl, $tag_txt);
    push @cmd_lines, ($run_bam2cov_n, $run_bam2cov_p, $run_cov2RNA_n, $run_cov2RNA_p,
    $run_catNP, $run_tag2pos, $run_sort2candi);
}

# Step 4. Merge RNAs from different libraries
push @cmd_lines, "\n$pTime \t \# Step 4. Merge RNAs from different libraries. \"  \n";
my @tag_files = ();
foreach my $lib (sort @reads_libs){
    my $tag_txt = catfile($Coverage_dir, "$name.$lib\_tag\.txt");
    push @tag_files, $tag_txt;
}
my $merge_tmp = catfile($Merge_dir, "$name\_merge.tmp");
my $merge_txt = catfile($Merge_dir, "$name\_merge.txt");
my $run_merge = join " ", ('perl', $merge4all_pl, @tag_files, '>', $merge_tmp);
my $run_txt2temp = join " ", ('perl', $sort2temp_pl, $merge_tmp);
my $run_temp2pos = join " ", ('perl', $getPosition_pl, '-i', "$merge_tmp\.temp", '-f', $gff, '-o', $merge_txt);
my $run_pos2candi = join " ", ('perl', $sort2candi_pl, $merge_txt);
push @cmd_lines, ($run_merge, $run_txt2temp, $run_temp2pos, $run_pos2candi);

# Step 5. Count reads
push @cmd_lines, "\n$pTime \t \# Step 5. Count reads of different libraries. \" ";
my $merge_gff = catfile($Count_dir, "$name\_merge.gff");
my $tmp_ptt = catfile($original_dir, 'tmp.ptt');
my $run_sort2ptt = join " ", ('perl', $sort2bed_pl, '-t sort2ptt', '-i', $merge_txt, '-o', $tmp_ptt);
my $run_ptt2gff = join " ", ('perl', $sort2bed_pl, '-t ptt2gff', '-i', $tmp_ptt, '-o', $merge_gff, '-g', $reference);
push @cmd_lines, ($run_sort2ptt, $run_ptt2gff);

foreach my $lib(sort @reads_libs){
    my $exp_out = catfile($Count_dir, "$name\_merge.lib$lib.txt");
    my $in_bam = catfile($Alignment_dir, "$name.$lib.s.bam");
    my $run_count = &Count_reads($merge_gff, $in_bam, $exp_out, 'multicov');
    push @cmd_lines, ($run_count);
}

# Statistic reads mapping (BAM files)
my $stat_map = catfile($original_dir, "$name\.mapping\.stat");
my $stat_log = catfile($original_dir, "$name\.mapping\.log");
my $run_stat = join " ", ('perl', $statMapping_pl, '-b', $Alignment_dir, '-o', $stat_map, '-l', $stat_log);

push @cmd_lines, "\n$pTime \t \# Finish. \"  \n";

print join "\n", @cmd_lines, "\n";


sub Count_reads{
    my ($in_gff, $in_bam, $out_txt, $type) = @_;
    $in_gff = abs_path($in_gff);
    $in_bam = abs_path($in_bam);
    $out_txt = abs_path($out_txt);
    my $htseq_count_cmd = join " ", ($htseq_count, '-f bam -o $out_txt', $in_bam, $in_gff, '>', $out_txt );
    my $multicov_cmd = join " ", ($bedtools, 'multicov -s -bams', $in_bam, '-bed', $in_gff, '>', $out_txt);
    if($type eq 'htseq-count'){
        return $htseq_count_cmd;
    }elsif($type eq 'multicov'){
        return $multicov_cmd;
    }else{
        return ;
    }
}

sub Check_file_existence{
    foreach my $f (@_){
        die "[$f]: File not exist." unless -e $f;
    }
}

sub check_dir{
    foreach my $i (@_){
        mkdir $i unless (-d $i);
    }
}

################################################################################
# Create PBS template
my $PBS_Header = "#PBS -S /bin/bash
PBS -q Debug
PBS -N Align_job.pbs
PBS -o $original_dir/job.out
PBS -e $original_dir/job.err

cd \$PBS_O_WORKDIR
#nohup  sh YOUR.sh > YOUR.sh.log &

# Write your Command here ... : Check the queue: Debug/General
# sh  \$PBS_O_WORKDIR/work.sh       > \$PBS_O_WORKDIR/work.sh.log

";

open H,"> wm_job.pbs" or die;
print H $PBS_Header;


__END__

=head1 NAME

c<Align.pl> - Create shell commands to mapping reads to reference.

=head1 SYNOPSIS

 # perl Align.pl -p 2 -r ../reads_clean -d ../database/H37Rv.fa -n H37Rv -x bowtie_clean > out.sh

=head1 OPTIONS

=over 8

=item B<-p> N, B<--Aligner>=N 

Choose the tools for the alignment.
1=soap, 2=bowtie, 3=bowtie2 [default: 1]

=item B<-r> str, B<--reads>=str

The dir of sequencing reads. 
The name of reads should consist of the following parts:
[Name][lib][part].[fa/fastq]

SE reads: [Name][01-02].[fa/fastq], eg: H37Rv01.fastq

PE reads: [Name][03-04]_[1-2].[fastq], eg: H37Rv03_1.fastq

=item B<-d> str, B<--reference>=str

The reference file in fasta format. single file.
Also need GFF file in the same directory of reference.fa.

eg: H37Rv.fa H37Rv.gff

=item B<-n> str, B<--name>=str

The name of this sample.
Using this word to search fa/fastq files in reads/ directory.

=item B<-x> str, B<--suffix>=str

The type of output dir. The output directory are: 
02.Alignmen_<suffix>, 03.Coverage_<suffix>, 04.Merge_tags_<suffix>
05.Count_reads_<suffix>, 06.RNAz_analysis_<suffix>

=item B<-h --help>

Print this help

=back

=head1 DESCRIPTION

This script will create the Shell commands to mapping reads to referenece.

=head1 AUTHOR

Wang Ming <wangmcas@gmail.com>

=cut
