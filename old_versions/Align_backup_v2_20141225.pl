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
&check_dir($Alignment_dir, $Coverage_dir, $Merge_dir, $Count_dir);

sub check_dir{
    foreach my $i (@_){
        mkdir $i unless (-d $i);
    }
}

## Check tools
my $soap        = which('soap');
my $bowtie      = which('bowtie');
my $bowtie2     = which('bowtie2');
my $samtools    = which('samtools');
my $bedtools    = which('bedtools');
my $htseq_count = which('htseq-count');
&Check_file_existence($soap, $bowtie, $bowtie2, $samtools, $bedtools, $htseq_count);

my $my_bin         = catdir($HOME, 'work/bin/get_sRNA');
my $get_sRNA_pl    = catfile($my_bin, 'get_sRNA.pl');
my $getPosition_pl = catfile($my_bin, 'getPosition_v4.pl');
my $merge4all_pl   = catfile($my_bin, 'merge4all.pl');
my $sort2candi_pl  = catfile($my_bin, 'sort2candi_v1.pl');
my $sort2temp_pl   = catfile($my_bin, 'sort2temp.pl');
my $sort2bed_pl    = catfile($my_bin, 'sort2bed.pl');
my $soapPE_merge_pl= catfile($my_bin, 'soapPE_merge.pl');
my $soap2sam_pl    = catfile($my_bin, 'soap2sam.pl');
&Check_file_existence($get_sRNA_pl, $getPosition_pl, $merge4all_pl, $sort2candi_pl,
                      $sort2temp_pl, $sort2bed_pl, $soapPE_merge_pl, $soap2sam_pl);

sub Check_file_existence{
    foreach my $f (@_){
        die "[$f]: File not exist." unless -e $f;
    }
}

## Check input Aligner
pod2usage(-verbose => 1, -message => "$0: -p | --Aligner, $Aligner not in [1-3]. \n") 
          unless ($Aligner == 1 || $Aligner == 2 || $Aligner == 3);

## Check reads
my $reads_dir = $reads;
my @list_reads = <$reads_dir/*.fa*>;
die "$0: -r | --reads, [$reads_dir], No fa/fastq files found. \n" if(@list_reads < 1);
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

# Step 1. Filter reads
push @cmd_lines, '#!/bin/bash';
push @cmd_lines, "$pTime \t \# Step 1. Filter reads \"  \n";

# Step 2. Create Alignment command lines
push @cmd_lines, "$pTime \t \# Step 2. Mapping reads and Find sRNAs. \"  \n";

$Alignment_dir = abs_path($Alignment_dir);
# Create soap commands
if($Aligner == 1){ # 1=soap
    my $soap_SE_para = '-v 2 -r 2 -p 2 -M 4';
    my $soap_PE_para = '-m 2 -r 2 -p 2 -x 500 -s 40 -l 32';
    # Check index file
    my $soap_index = $reference.".index";
    my @check_index = <$soap_index.*>;
    my $soap_built = which('2bwt-builder');
    system "$soap_build $reference" unless (@check_index < 1);

    foreach my $lib (@reads_libs){
        push @cmd_lines, "$pTime \t \#\#\# Mapping reads: $name $lib \"  \n ";
        my $soap_out = catfile($Alignment_dir, "$name.$lib.soap");
        my $soap_PE = catfile($Alignment_dir, "$name.$lib.soapPE");
        my $soap_PESingle = catfile($Alignment_dir, "$name.$lib.soapSingle");
        my $soap_un = catfile($Alignment_dir, "$name.$lib.unsoap");
# Run soap
        # Soap for SE
        if($lib eq "01" || $lib eq "02"){
            my $SE_reads_name = $reads_info{$lib}[0];
            my $SE_reads = abs_path(catfile($reads_dir, $SE_reads_name));
            my $run_soap = join " ", ($Aligner $soap_SE_para '-D', $soap_index, 
                    '-a', $SE_reads, '-o', $soap_out, '-u', $soap_un);
            push @cmd_lines, $run_soap;
        # Soap for PE
        }elsif($lib >= "03"){
            my @PE_reads = @reads_info{$lib};
            my $PE_reads_1 = abs_path(catfile($reads_dir, $PE_reads[0]));
            my $PE_reads_2 = abs_path(catfile($reads_dir, $PE_reads[1]));
            my $run_soap = join " ", ($Aligner $soap_PE_para, '-D', $soap_index, 
                    '-a', $PE_reads_1, '-b', $PE_reads_2, '-o' $soap_out, 
                    '-2', $soap_PESingle, '-u', $soap_un);
            my $soap_PE2SE_fa  = catfile($Alignment_dir, "$name.$lib.fa");
            my $soap_PE2SE_del = catfile($Alignment_dir, "$name.$lib.del");
            my $run_soapPEmerge = join " ", ('perl', $soapPE_merge_pl, '-a', $reference, 
                    '-i', $soap_PE, '-o', $soap_out, '-c', $soap_PE2SE_fa, '-d',
                    $soap_PE2SE_del );
            push @cmd_lines, ($run_soap, $run_soapPEmerge);
        }
# Run soap to bam
        my $sam_prefix = catfile($Alignment_dir, "$name.$lib");
        my $sam_file = $sam_prefix.'.sam';
        my $bam_file = $bam_prefix.'.bam';
        my $bam_s_file = $bam_prefix.'.s.bam';
        my $run_soap2sam = join " ", ('perl', $soap2sam, $soap_out, '>', $sam_file);
        my $run_sam2bam  = join " ", ($samtools, 'view -bt', "$reference\.fai", $sam_file, 
                '-o', $bam_file);
        my $run_bam2sort = join " ", ($samtools, 'sort', $bam_file, "$sam_prefix\.s");
        my $run_bam2index = join " ", ($samtools, 'index', $bam_s_file);
        push @cmd_lines, ($run_soap2sam, $run_sam2bam, $run_bam2sort, $run_bam2index);

# Run bam to coverage
        push @cmd_lines, "$pTime \t \#\#\# Find sRNA. $name $lib\"  \n ";
        my $Coverage_n = catfile($Coverage_dir, "$name.$lib.coverage.n");
        my $Coverage_p = catfile($Coverage_dir, "$name.$lib.coverage.p");
        my $tag_tmp_n  = catfile($Coverage_dir, "$name.$lib\_tag\.tmp\.n");
        my $tag_tmp_p  = catfile($Coverage_dir, "$name.$lib\_tag\.tmp\.p");
        my $tag_tmp    = catfile($Coverage_dir, "$name.$lib\_tag\.tmp");
        my $tag_txt    = catfile($Covearge_dir, "$name.$lib\_tag\.txt");
        my $run_bam2cov_n = join " ", ($bedtools, 'genomecov', '-ibam', $bam_s_file, '-d -g', 
                "$reference\.fai", '-strand - >', $Coverage_n);
        my $run_bam2cov_p = join " ", ($bedtools, 'genomecov', '-ibam', $bam_s_file, '-d -g', 
                "$reference\.fai", '-strand + >', $Coverage_p);
        my $run_cov2RNA_n = join " ", ('perl', $get_sRNA_pl, '-i', $Coverage_n, 
                '-c 0.5 -s - -o', $tag_tmp_n);
        my $run_cov2RNA_p = join " ", ('perl', $get_sRNA_pl, '-i', $Coverage_p,
                '-c 0.5 -s - -o', $tag_tmp_p);
        my $run_catNP = join " ", ('cat', $tag_tmp_n, $tag_tmp_p, '>', $tag_tmp);
        my $run_tag2pos = join " ", ('perl', $getPosition_pl, '-i', $tag_tmp, '-f', $gff, 
                '-o', $tag_txt);
        my $run_sort2candi = join " ", ('perl', $sort2candi_pl, $tag_txt);
        push @cmd_lines, ($run_bam2cov_n, $run_bam2cov_p, $run_cov2RNA_n, $run_cov2RNA_p,
                $run_catNP, $run_tag2pos, $run_sort2candi);
    }
# Create bowtie commands
}elsif($Aligner == 2){ # bowtie
    my $bowtie_SE_para = '-v 3 --best --strata -S -m 100 -X 300 --chunkmbs 256 -p 4';
    my $bowtie_PE_para = '-v 3 --best --strata -S -m 100 -I 60 -X 300 --chunkmbs 256 -p 4';
    # Check index file
    my $reference_dir = dirname($reference);
    my $bowtie_index = catfile($reference_dir, $name);
    my @check_index  = <$bowtie_index\.*\.ebwt>;
    my $bowtie_build = which('bowtie-build'); &Check_file_existence($bowtie_build);
    system"$bowtie_build -q -f $reference $bowtie_index" if(@check_index < 1);
    
    foreach my $lib (@reads_libs){
        push @cmd_lines, "$pTime, \t \#\#\# Mappping reads: $name $lib \" \n ";
        my $sam_prefix = catfile($Alignment_dir, "$name.$lib");
        my $sam_file   = catfile($Alignment_dir, "$name.$lib.sam");
        my $bam_file   = catfile($Alignment_dir, "$name.$lib.bam");
        my $bam_s_file = catfile($Alignment_dir, "$name.$lib.s.bam");
# Run bowtie (sam output)
        my $run_bowtie = '';
        my $read_format = '-q';
        if($lib eq '01' || $lib eq '02'){
            my $SE_reads_name = $reads_info{$lib}[0];
            my $SE_reads = abs_path(catfile($reads_dir, $SE_reads_name));
            $read_format = '-f' if($SE_reads =~ /\.fa$/);
            $run_bowtie = join " ", ($bowtie, $bowtie_SE_para, $bowtie_index, $read_format,
                    $SE_reads, $sam_file);
        }elsif($lib >= '03'){
            my @PE_reads = @reads_info{$lib};
            my $PE_reads_1 = abs_path(catfile($reads_dir, $PE_reads[0]));
            my $PE_reads_2 = abs_path(catfile($reads_dir, $PE_reads[1]));
            $read_format = '-f' if($PE_reads_1 =~ /\.fa$/);
            $run_bowtie = join " ", ($bowtie, $bowtie_PE_para, $bowtie_index, 
                    '-1', $PE_reads_1, '-2', $PE_reads_2, $sam_file);        
        }else{
            die "Check the type of reads: $lib.";
        }
        my $run_sam2bam = join " ", ($samtools, 'view -bt', "$reference\.fai", $sam_file,
                '-o', $bam_file);
        my $run_bam2sort = join " ", ($samtools, 'sort', $bam_file, "$sam_prefix.s");
        my $run_bam2index = join " ", ($samtools, 'index', $bam_s_file);
        push @cmd_lines, ($run_bowtie, $run_sam2bam, $run_bam2sort, $run_bam2index);
# Run bam to coverage, Find RNA
        push @cmd_lines, "$pTime \t \#\#\# Find sRNA. $name $lib\"  \n ";
        my $Coverage_n = catfile($Coverage_dir, "$name.$lib.coverage.n");
        my $Coverage_p = catfile($Coverage_dir, "$name.$lib.coverage.p");
        my $tag_tmp_n  = catfile($Coverage_dir, "$name.$lib\_tag\.tmp\.n");
        my $tag_tmp_p  = catfile($Coverage_dir, "$name.$lib\_tag\.tmp\.p");
        my $tag_tmp    = catfile($Coverage_dir, "$name.$lib\_tag\.tmp");
        my $tag_txt    = catfile($Covearge_dir, "$name.$lib\_tag\.txt");
        my $run_bam2cov_n = join " ", ($bedtools, 'genomecov', '-ibam', $bam_s_file, '-d -g', 
                "$reference\.fai", '-strand - >', $Coverage_n);
        my $run_bam2cov_p = join " ", ($bedtools, 'genomecov', '-ibam', $bam_s_file, '-d -g', 
                "$reference\.fai", '-strand + >', $Coverage_p);
        my $run_cov2RNA_n = join " ", ('perl', $get_sRNA_pl, '-i', $Coverage_n, 
                '-c 0.5 -s - -o', $tag_tmp_n);
        my $run_cov2RNA_p = join " ", ('perl', $get_sRNA_pl, '-i', $Coverage_p,
                '-c 0.5 -s - -o', $tag_tmp_p);
        my $run_catNP = join " ", ('cat', $tag_tmp_n, $tag_tmp_p, '>', $tag_tmp);
        my $run_tag2pos = join " ", ('perl', $getPosition_pl, '-i', $tag_tmp, '-f', $gff, 
                '-o', $tag_txt);
        my $run_sort2candi = join " ", ('perl', $sort2candi_pl, $tag_txt);
        push @cmd_lines, ($run_bam2cov_n, $run_bam2cov_p, $run_cov2RNA_n, $run_cov2RNA_p,
                $run_catNP, $run_tag2pos, $run_sort2candi);    
    }
# Create bowtie2 commands
}elsif($Aligner == 3){ # bowtie2
    my $bowtie2_SE_para = '';
    my $bowtie2_PE_para = '';
    # Check index file
    my $reference_dir = dirname($reference);
    my $bowtie2_index = catfile($reference_dir, $name);
    my @check_index = <$bowtie2_index\..*\.ebwt>;
    my $bowtie2_build = which('bowtie2-build'); &Check_file_existence($bowtie2_build);
    system"$bowtie2_build -q -f $reference $bowtie2_index";

    foreach my $lib (@reads_libs){
        push @cmd_lines, "$pTime, \t \#\#\# Mappping reads: $name $lib \" \n ";
        my $sam_prefix = catfile($Alignment_dir, "$name.$lib");
        my $sam_file   = catfile($Alignment_dir, "$name.$lib.sam");
        my $bam_file   = catfile($Alignment_dir, "$name.$lib.bam");
        my $bam_s_file = catfile($Alignment_dir, "$name.$lib.s.bam");
# Run bowtie2 (sam output)
        my $run_bowtie2 = '';
        my $read_format = '-q';
        if($lib eq '01' || $lib eq '02'){
            my $SE_reads_name = $reads_info{$lib}[0];
            my $SE_reads = abs_path(catfile($reads_dir, $SE_reads_name));
            $read_format = '-f' if($SE_reads =~ /\.fa$/);
            $run_bowtie2 = join " ", ($bowtie2, $bowtie2_SE_para, '-x', $bowtie2_index, 
                    $read_format, $SE_reads, '-S', $sam_file);
        }elsif($lib >= '03'){
            my @PE_reads = @reads_info{$lib};
            my $PE_reads_1 = abs_path(catfile($reads_dir, $PE_reads[0]));
            my $PE_reads_2 = abs_path(catfile($reads_dir, $PE_reads[1]));
            $read_format = '-f' if($PE_reads_1 =~ /\.fa$/);
            $run_bowtie2 = join " ", ($bowtie2, $bowtie2_PE_para, '-x', $bowtie2_index,
                    '-1', $PE_reads_1, '-2', $PE_reads_2, '-S', $sam_file);        
        }else{
            die "Check the type of reads: $lib.";
        }
        my $run_sam2bam = join " ", ($samtools, 'view -bt', "$reference\.fai", $sam_file,
                '-o', $bam_file);
        my $run_bam2sort = join " ", ($samtools, 'sort', $bam_file, "$sam_prefix.s");
        my $run_bam2index = join " ", ($samtools, 'index', $bam_s_file);
        push @cmd_lines, ($run_bowtie, $run_sam2bam, $run_bam2sort, $run_bam2index);
# Run bam to coverage, Find RNA
        push @cmd_lines, "$pTime \t \#\#\# Find sRNA. $name $lib\"  \n ";
        my $Coverage_n = catfile($Coverage_dir, "$name.$lib.coverage.n");
        my $Coverage_p = catfile($Coverage_dir, "$name.$lib.coverage.p");
        my $tag_tmp_n  = catfile($Coverage_dir, "$name.$lib\_tag\.tmp\.n");
        my $tag_tmp_p  = catfile($Coverage_dir, "$name.$lib\_tag\.tmp\.p");
        my $tag_tmp    = catfile($Coverage_dir, "$name.$lib\_tag\.tmp");
        my $tag_txt    = catfile($Covearge_dir, "$name.$lib\_tag\.txt");
        my $run_bam2cov_n = join " ", ($bedtools, 'genomecov', '-ibam', $bam_s_file, '-d -g', 
                "$reference\.fai", '-strand - >', $Coverage_n);
        my $run_bam2cov_p = join " ", ($bedtools, 'genomecov', '-ibam', $bam_s_file, '-d -g', 
                "$reference\.fai", '-strand + >', $Coverage_p);
        my $run_cov2RNA_n = join " ", ('perl', $get_sRNA_pl, '-i', $Coverage_n, 
                '-c 0.5 -s - -o', $tag_tmp_n);
        my $run_cov2RNA_p = join " ", ('perl', $get_sRNA_pl, '-i', $Coverage_p,
                '-c 0.5 -s - -o', $tag_tmp_p);
        my $run_catNP = join " ", ('cat', $tag_tmp_n, $tag_tmp_p, '>', $tag_tmp);
        my $run_tag2pos = join " ", ('perl', $getPosition_pl, '-i', $tag_tmp, '-f', $gff, 
                '-o', $tag_txt);
        my $run_sort2candi = join " ", ('perl', $sort2candi_pl, $tag_txt);
        push @cmd_lines, ($run_bam2cov_n, $run_bam2cov_p, $run_cov2RNA_n, $run_cov2RNA_p,
                $run_catNP, $run_tag2pos, $run_sort2candi);
    }
# Other type of Aligner
}else{
    die "$0: -p | --Aligner [$Aligner] Not in [soap, bowtie, bowtie2] \n";
}

################################################################################
## SOAP commands
#    my $soap_cmd_SE = " -v 2 -r 2 -p 2 -M 4 ";
#    my $soap_cmd_PE = " -m 2 -r 2 -p 2 -x 500 -s 40 -l 32 ";
#    foreach my $lib (@libs){
###############################
## A. Create SOAP command lines
#        push @Align_cmd_lines, "\n$pTime \t \#\#\# Mapping: $opts{n}$lib \" ";
#        my @Reads2bam_cmd = ();
#        my $lib_out_prefix = $Alignment_path."/".$opts{n}.$lib;
#        if($lib eq "01" || $lib eq "02"){
#            my @reads = <$reads_abs_path\/$opts{n}$lib\.fa*>;  # search SE fa/fastq
#            die "Find more $reads_abs_path\/$opts{n}\.$lib\. files " if(@reads > 1); ###
#            my $lib_reads = shift(@reads);
##            my $lib_reads = $reads_abs_path."/".$opts{n}.$lib.".fa";  ##
##            my $lib_out_prefix = $Alignment_path."/".$opts{n}.$lib;
#            my $R_soap = $lib_out_prefix.".soap";
#            my $R_unsoap = $lib_out_prefix.".unsoap";        
#            $Reads2bam_cmd[0] = join" ",($tools{Aligner}, $soap_cmd_SE, "-D", $Aligner_index, "-a", $lib_reads, "-o", $R_soap, "-u", $R_unsoap );
#            &Check_file_exists($tools{Aligner}, $lib_reads);
###############################
## 1). SE soap -> sam
#            my $soapSE_sam = $lib_out_prefix.".sam";
#            my $soapSE_bam = $lib_out_prefix.".bam";
#            my $soapSE_bam_sorted = $lib_out_prefix.".s.bam";
#            $Reads2bam_cmd[1] = "$tools{perl} $tools{soap2sam_pl}  $R_soap  > $soapSE_sam";
#            $Reads2bam_cmd[2] = "$tools{samtools} view -bt $Ref_genome\.fai $soapSE_sam -o $soapSE_bam";
#            $Reads2bam_cmd[3] = "$tools{samtools} sort $soapSE_bam $lib_out_prefix\.s";
#            $Reads2bam_cmd[4] = "$tools{samtools} index $soapSE_bam_sorted";
#            push @Align_cmd_lines, @Reads2bam_cmd;
##        }elsif($lib eq "03" || $lib eq "04") {
#        }elsif($lib > "02"){
#            my $lib_reads_1 = $reads_abs_path."/".$opts{n}.$lib."_1.fastq";
#            my $lib_reads_2 = $reads_abs_path."/".$opts{n}.$lib."_2.fastq";
##            my $lib_out_prefix = $Alignment_path."/".$opts{n}.$lib;
#            my $R_soap = $lib_out_prefix.".PEsoap";
#            my $R_unsoap = $lib_out_prefix.".PEunsoap";
#            my $R_Singlesoap = $lib_out_prefix.".Singlesoap";
#            $Reads2bam_cmd[0] = join" ",($tools{Aligner}, $soap_cmd_PE, "-D", $Aligner_index, "-a", $lib_reads_1, "-b", $lib_reads_2, "-o", $R_soap, "-2", $R_Singlesoap, "-u", $R_unsoap);
#            &Check_file_exists($tools{Aligner}, $lib_reads_1, $lib_reads_2);
###############################
## 2). PE soap -> sam
#            my $soapPE_merge = $lib_out_prefix.".soap";
#            my $soapPE_merge_fa = $lib_out_prefix.".fa";
#            my $soapPE_merge_del = $lib_out_prefix.".del";
#            my $soapPE_merge_sam = $lib_out_prefix.".sam";
#            my $soapPE_merge_bam = $lib_out_prefix.".bam";
#            my $soapPE_merge_bam_sorted = $lib_out_prefix.".s.bam";
#            $Reads2bam_cmd[1] = "$tools{perl} $tools{soapPE_merge_pl} -a $Ref_genome -i $R_soap -o $soapPE_merge -c $soapPE_merge_fa -d $soapPE_merge_del";
#            $Reads2bam_cmd[2] = "$tools{perl} $tools{soap2sam_pl} $soapPE_merge  > $soapPE_merge_sam";
#            $Reads2bam_cmd[3] = "$tools{samtools} view -bt $Ref_genome\.fai $soapPE_merge_sam -o $soapPE_merge_bam";
#            $Reads2bam_cmd[4] = "$tools{samtools} sort $soapPE_merge_bam  $lib_out_prefix\.s ";
#            $Reads2bam_cmd[5] = "$tools{samtools} index $soapPE_merge_bam_sorted";
#            push @Align_cmd_lines, @Reads2bam_cmd;
#        }else{
#            die "Check the value of @libs, should be 01-04 ";
#        }
##        push @Align_cmd_lines,"\n";
#    }
#} elsif($opts{p} eq "bowtie") {
################################################################################
## B. Create Bowtie commands
#    my $bowtie_cmd_SE = " -v 3 --best --strata -S -m 100 -X 300 --chunkmbs 256 -p 4 ";
#    my $bowtie_cmd_PE = " -v 3 --best --strata -S -m 100 -I 60 -X 300 --chunkmbs 256 -p 4 ";
#    # SE  bowtie -para index reads
#    foreach my $lib (@libs){
#        my $lib_out_prefix = $Alignment_path."/".$opts{n}.$lib;
#        my $Out_sam = $lib_out_prefix.".sam";
#        push @Align_cmd_lines,"\n\n$pTime \t \#\#\# Mapping: $opts{n}$lib\" ";
#        if($lib eq "01" || $lib eq "02"){
###############################
## 1). SE => sam
#            my @reads = <$reads_abs_path\/$opts{n}$lib\.fa*>;
#            die "More than one files found: $reads_abs_path\/$opts{n}\.$lib\.fa* " if(@reads > 1);
#            my $lib_reads = shift(@reads);
#            my $reads_para = "-f";
#            $reads_para = "-q" if($lib_reads =~ /\.fastq$/);
#            my $SE_cmd = join " ", ($tools{Aligner}, $bowtie_cmd_SE, $Aligner_index, $reads_para, $lib_reads, $Out_sam);
##            my $lib_reads = $reads_abs_path."/".$opts{n}.$lib.".fa";
##            my $SE_cmd = join" ", ($tools{Aligner}, $bowtie_cmd_SE, $Aligner_index, "-f", $lib_reads, $Out_sam);
#            &Check_file_exists($tools{Aligner}, $lib_reads);
#            push @Align_cmd_lines, $SE_cmd;
##        }elsif($lib eq "03" || $lib eq "04"){
#        }elsif($lib > "02"){
###############################
## 2). PE => sam
#            my $lib_reads_1 = $reads_abs_path."/".$opts{n}.$lib."_1.fastq";
#            my $lib_reads_2 = $reads_abs_path."/".$opts{n}.$lib."_2.fastq";
#            my $PE_cmd = join" ", ($tools{Aligner}, $bowtie_cmd_PE, $Aligner_index, "-1", $lib_reads_1, "-2", $lib_reads_2, $Out_sam);
#            &Check_file_exists($tools{Aligner}, $lib_reads_1, $lib_reads_2);
#            push @Align_cmd_lines, $PE_cmd;
#        }else{
#            die "Check the value of @libs, should be 01-04 ";
#        }
#        my $Out_bam = $lib_out_prefix.".bam";
#        my $Out_bam_sorted = $lib_out_prefix.".s.bam";
#        my $sam_cmd_1 = "$tools{samtools} view -bt $Ref_genome\.fai $Out_sam -o $Out_bam";
#        my $sam_cmd_2 = "$tools{samtools} sort $Out_bam $lib_out_prefix\.s";
#        my $sam_cmd_3 = "$tools{samtools} index $Out_bam_sorted";
#        push @Align_cmd_lines, ($sam_cmd_1, $sam_cmd_2, $sam_cmd_3);
#        push @Align_cmd_lines,"\n";
#    }
#}elsif($opts{p} eq "bowtie2"){
#################################################################################
## C. Create Bowtie2 commands
#    my $bowtie2_cmd_SE = "-p 4";
#    my $bowtie2_cmd_PE = "-I 50 -X 200 -p 4";
#    foreach my $lib (@libs){
#        my $lib_out_prefix = $Alignment_path."/".$opts{n}.$lib;
#        my $Out_sam = $lib_out_prefix.".sam";
#        push @Align_cmd_lines, "\n\n$pTime \t \#\#\# Mapping: $opts{n}$lib\" ";
#        if($lib eq "01" || $lib eq "02"){
###############################
## 1). SE => sam
##            my $lib_reads = $reads_abs_path."/".$opts{n}.$lib.".fa";
##            my $SE_cmd = join" ", ($tools{Aligner}, $bowtie2_cmd_SE, "-x", $Aligner_index, "-f", $lib_reads, "-S", $Out_sam);
#            my @reads = <$reads_abs_path\/$opts{n}$lib\.fa*>;
#            die "More than one files found: $reads_abs_path\/$opts{n}\.$lib\.fa* " if(@reads > 1);
#            my $lib_reads = shift(@reads);
#            my $reads_para = "-f";
#            $reads_para = "-q" if($lib_reads =~ /\.fastq$/);
#            my $SE_cmd = join " ", ($tools{Aligner}, $bowtie2_cmd_SE, "-x", $Aligner_index, $reads_para, $lib_reads, "-S", $Out_sam);
#            &Check_file_exists($tools{Aligner}, $lib_reads);
#            push @Align_cmd_lines, $SE_cmd;
##        }elsif($lib eq "03" || $lib eq "04"){
#        }elsif($lib > "02"){
###############################
## 2). PE => sam            
#            my $lib_reads_1 = $reads_abs_path."/".$opts{n}.$lib."_1.fastq";
#            my $lib_reads_2 = $reads_abs_path."/".$opts{n}.$lib."_2.fastq";
#            my $PE_cmd = join" ", ($tools{Aligner}, $bowtie2_cmd_PE, "-x", $Aligner_index, "-1", $lib_reads_1, "-2", $lib_reads_2, "-S", $Out_sam);
#            &Check_file_exists($tools{Aligner}, $lib_reads_1, $lib_reads_2);
#            push @Align_cmd_lines, $PE_cmd;
#        }else{
#            die "Check the value of @libs, should be 01-04 ";
#        }
#        my $Out_bam = $lib_out_prefix.".bam";
#        my $Out_bam_sorted = $lib_out_prefix.".s.bam";
#        my $sam_cmd_1 = "$tools{samtools} view -bt $Ref_genome\.fai $Out_sam -o $Out_bam";
#        my $sam_cmd_2 = "$tools{samtools} sort $Out_bam $lib_out_prefix\.s";
#        my $sam_cmd_3 = "$tools{samtools} index $Out_bam_sorted";
#        push @Align_cmd_lines, ($sam_cmd_1, $sam_cmd_2, $sam_cmd_3);
#        push @Align_cmd_lines,"\n";
#    }
#} else {
#    die "Check \$opts{p}, should be soap or bowtie";
#}
##print join"\n",@Align_cmd_lines,"\n";
#
#sub Check_file_exists{
#    foreach(@_){
#        die "File Not Found:\n$_" unless -e $_;
#    }
#}
#
#################################################################################
## Step 3. BAM => coverage => sRNA
## Create Coverage dir
#my $Coverage_dir = $opts{C};
#mkdir $Coverage_dir unless -d $Coverage_dir;
#my $Coverage_path = abs_path($Coverage_dir);
#
## BAM to Coverage
#push @Align_cmd_lines, "\n$pTime \t \# Step 3. Transfer BAM to Coverage file \" ";
#foreach my $lib (@libs) {
#    my $bam_in = $work_dir."/".$opts{A}."/".$opts{n}.$lib.".s.bam";
##    &Check_file_exists($bam_in);
#    my $bam_in_abs_path = abs_path($bam_in);
#    my $bam_name = basename($bam_in);
#    my $Coverage_prefix = $opts{n}.$lib;
#    my @bam2cov_cmd = ();
#    $bam2cov_cmd[0] = "\n\# Sample: $opts{n}$lib ";
#    $bam2cov_cmd[1] = "$tools{genomeCoverageBed} -ibam $bam_in_abs_path -d -g $opts{n}\.fai -strand - > $Coverage_path\/$Coverage_prefix\.coverage\.n ";
#    $bam2cov_cmd[2] = "$tools{genomeCoverageBed} -ibam $bam_in_abs_path -d -g $opts{n}\.fai -strand + > $Coverage_path\/$Coverage_prefix\.coverage\.p ";
#    $bam2cov_cmd[3] = "$tools{perl} $tools{get_sRNA_pl} -i $Coverage_path\/$Coverage_prefix\.coverage\.n -c 0.5 -s - -o $Coverage_path\/$Coverage_prefix\_tag.tmp\.n ";
#    $bam2cov_cmd[4] = "$tools{perl} $tools{get_sRNA_pl} -i $Coverage_path\/$Coverage_prefix\.coverage\.p -c 0.5 -s + -o $Coverage_path\/$Coverage_prefix\_tag.tmp\.p ";
#    $bam2cov_cmd[5] = "cat $Coverage_path\/$Coverage_prefix\_tag.tmp\.n $Coverage_path\/$Coverage_prefix\_tag.tmp\.p > $Coverage_path\/$Coverage_prefix\_tag.tmp ";
#    $bam2cov_cmd[6] = "$tools{perl} $tools{getPosition_pl} -i $Coverage_path\/$Coverage_prefix\_tag.tmp -f $opts{f} -o $Coverage_path\/$Coverage_prefix\_tag\.txt ";
#    $bam2cov_cmd[7] = "$tools{perl} $tools{sort2candi_pl}  $Coverage_path\/$Coverage_prefix\_tag\.txt ";
#    push @Align_cmd_lines, @bam2cov_cmd;
#}
#
#################################################################################
## Step 4. Merge sRNAs from all libraries
#mkdir $opts{M} unless -d $opts{M}; # 04.Merge_tags
#my $merge_path = abs_path($opts{M});
#
#push @Align_cmd_lines, "\n$pTime \t \# Step 4. Merge sRNA tags from all libraries\" \n";
#my @tag_files;
#foreach my $lib (@libs){
#    my $t_file = "$Coverage_path\/$opts{n}$lib\_tag\.txt";
#    push @tag_files, $t_file;
#}
#my $parameter_file = join" ",(@tag_files);
#my $merge_cmd = "$tools{perl} $tools{merge4all_pl}  $parameter_file  > $merge_path\/$opts{n}\_merge.tmp ";
#my $merge2temp = "$tools{perl} $tools{sort2temp_pl} $merge_path\/$opts{n}\_merge.tmp";
#my $merge2temp2pos = "$tools{perl} $tools{getPosition_pl}  -i $merge_path\/$opts{n}\_merge.tmp.temp -f $opts{f} -o $merge_path\/$opts{n}\_merge.txt";
#my $merge2sort = "$tools{perl} $tools{sort2candi_pl}  $merge_path\/$opts{n}\_merge.txt";
#push @Align_cmd_lines, ($merge_cmd, $merge2temp, $merge2temp2pos, $merge2sort);
#
#push @Align_cmd_lines, "\n$pTime \t \# Finsh. Write sRNAs to 04.Merge/$opts{n}\_merge.txt \" ";
#
##print join"\n",@Align_cmd_lines,"\n";
#
#################################################################################
## Step 5. Count reads on each sRNAs
## 1. htset-count -f bam  --stranded=yes -q Alignment.s.bam in.gff > out.exp
## 2. bedtools multicov -s -bams align1.bam align2.bam ... -bed in.gff >out.exp
#push @Align_cmd_lines, "$pTime \t \# Step 5. Count reads on each sRNAs\" ";
#
#mkdir $opts{E} unless -d $opts{E};
#my $count_path = abs_path($opts{E});
#
#my @Step5cmd = ();
#my $sRNA_file = abs_path("$merge_path\/$opts{n}_merge\.txt");
#my $sRNA_gff  = basename($sRNA_file);
#$sRNA_gff =~ s/\.txt/.gff/;
#$sRNA_gff = $count_path."/".$sRNA_gff;
#my $trans_sort2gff = "$tools{perl} $tools{sort2bed_pl} -t sort2gff -i $sRNA_file -o $sRNA_gff  -g $Ref_genome ";
#my $extra_1 = "head $sRNA_gff > tmp.gff";
#my $extra_2 = "mv -f tmp.gff $sRNA_gff";
#push @Align_cmd_lines, $trans_sort2gff;
#push @Align_cmd_lines,"\n";
##push @Step5cmd, ($trans_sort2gff, $extra_1, $extra_2, "\n");
#
##&Check_file_exists($sRNA_file); ## Update Nov-17
#
#foreach my $lib (@libs){
#    my $sorted_bam = "$Alignment_path\/$opts{n}$lib\.s\.bam";
##    warn "BAM file Not Found:\n$sorted_bam" unless -e $sorted_bam; ## Update Nov-17
#    my $lib_cmd = &CountReads($sRNA_file, $sorted_bam, $lib);
#    push @Align_cmd_lines, "$pTime \t \# Count lib: $lib \" ";
#    push @Align_cmd_lines, $lib_cmd;
##    push @Step5cmd, ("$pTime  \# Run lib: $lib \" ", $lib_cmd);
#}
#
##print join"\n",@Step5cmd;
#print join"\n", @Align_cmd_lines, "\n";
#
#sub CountReads{
#    my ($in_txt, $in_bam, $lib_id) = @_;
#    $in_txt = abs_path($in_txt);
#    $in_bam = abs_path($in_bam);
## Prepare output file names    
#    my $in_txt_exp = basename($in_txt);
#    $in_txt_exp =~ s/txt$/lib$lib_id\.txt/;
#    $in_txt_exp = $count_path."/".$in_txt_exp;
## Create commands
#    my $multicov_cmd = "$tools{bedtools} multicov -s -bams $in_bam -bed $sRNA_gff > $in_txt_exp ";
#    return join"\n",($multicov_cmd, "\n");
##    my $htseq_cmd   = "$tools{htseq_count}  -f bam -q -o out.sam  $in_bam  $sRNA_gff  > $in_txt_exp";    
##    return join"\n",($htseq_cmd,"\n");  
#}
#
#################################################################################
## Create PBS template
#my $PBS_Header = "#PBS -S /bin/bash
##PBS -q Debug
##PBS -o $work_dir/pbs.out
##PBS -e $work_dir/pbs.err
##
#cd \$PBS_O_WORKDIR
##nohup  sh YOUR.sh > YOUR.sh.log &
#
## Write your Command here ... : Check the queue: Debug/General
## sh  \$PBS_O_WORKDIR/work.sh       > \$PBS_O_WORKDIR/work.sh.log
##
#";
#
#open H,"> wm_job.pbs" or die;
#print H $PBS_Header;



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
