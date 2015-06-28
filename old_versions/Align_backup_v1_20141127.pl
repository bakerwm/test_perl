#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use Cwd 'abs_path';
use File::Basename  qw/basename dirname/;
# History:2014-11-23  updated: line-199, 245, 281; elsif($lib eq "03" || $lib eq "04")  => elsif($lib > "02")
sub usage{
    print << "EOF";
Usage: Create a shell script to Align short reads to reference genome 
using Bowtie/SOAP.

    -p :    Select the aligner, Bowtie/Bowtie2/SOAP, default: Bowtie
    -r :    The Path of short reads, [eg: /home/wangming/work/reads]
    -d :    The file of reference, [eg: *.fa, *.fna]
    -f :    The GFF of reference, [eg: *.gff]
    -n :    The name of the sample, [eg: H37Rv]
    -A :    The dir to store Alignment results, default: [02.Alignment]
    -C :    The dir to store Coverage files, default: [03.Coverage]
    -M :    The dir to store Merge sRNA file, default: [04.Merge_tags]
    -E :    The dir to store Count reads file, default: [05.Count_reads]
    -h :	Show this help

Example:
perl Align.pl -p soap -r reads -d NC_000962.fna -f NC_000962.gff -n H37Rv -A 02.Alignment -C 03.Coverage

Note:
1. Create directories: 01.Filter_reads, 02.Alignment, 03.Coverage, 04.Merge_tags
EOF
exit(1);
}

################################################################################
# Step 00. input parameters
my %opts = ();
getopts("p:r:d:f:n:A:C:M:E:h:", \%opts);
&usage if(! defined ($opts{r}) || ! defined ($opts{d}) || ! defined ($opts{f}) || defined ($opts{h}));
$opts{p} = "bowtie" unless (defined $opts{p});
$opts{n} = "temp" unless(defined $opts{n});
$opts{A} = "02.Alignment" unless (defined $opts{A});
$opts{C} = "03.Coverage" unless (defined $opts{C});
$opts{M} = "04.Merge_tags" unless (defined $opts{M});
$opts{E} = "05.Count_reads" unless (defined $opts{E});
die "-p: input should be bowtie/soap" unless($opts{p} =~ /^bowtie/ || $opts{p} eq "soap");

################################################################################
# Step 0. Prepare tools
# bowtie, bowtie-build, soap, 2bwt-builder,
# samtools, BEDTools, get_sRNA.pl, getPosition.pl sort2cand.pl 
my %tools = ();
#$tools{Aligner}           = `which $opts{p}`; chomp($tools{Aligner});
#$tools{Aligner_build}     = $tools{Aligner}."-build";
#if($tools{Aligner} eq "soap"){$tools{Aligner_build} = `which 2bwt-builder`; chomp($tools{Aligner_build})};
#$tools{samtools}          = `which samtools`; chomp($tools{samtools});
#$tools{genomeCoverageBed} = `which genomeCoverageBed`; chomp($tools{genomeCoverage});
if($opts{p} eq "soap"){
    $tools{Aligner}       = "/home/wangming/software/soap2.20release/soap";
    $tools{Aligner_build} = "/home/wangming/software/soap2.20release/2bwt-builder";
}elsif($opts{p} eq "bowtie"){
    $tools{Aligner}       = "/home/wangming/software/bowtie-1.1.1/bowtie";
    $tools{Aligner_build} = "/home/wangming/software/bowtie-1.1.1/bowtie-build"; 
}elsif($opts{p} eq "bowtie2"){
    $tools{Aligner}       = "/home/wangming/software/bowtie2-2.2.4/bowtie2";
    $tools{Aligner_build} = "/home/wangming/software/bowtie2-2.2.4/bowtie2-build";
}else{
    die "Choose a aligner: soap/bowtie/bowtie2";
}
$tools{samtools}          = "/home/wangming/software/samtools-1.1/samtools";
$tools{genomeCoverageBed} = "/home/wangming/software/bedtools2-2.21.0/bin/genomeCoverageBed";
$tools{get_sRNA_pl}       = "/home/wangming/work/bin/get_sRNA/get_sRNA.pl";
$tools{getPosition_pl}    = "/home/wangming/work/bin/get_sRNA/getPosition_v4.pl";
$tools{merge4all_pl}      = "/home/wangming/work/bin/get_sRNA/merge4all.pl";
$tools{sort2candi_pl}     = "/home/wangming/work/bin/get_sRNA/sort2candi_v1.pl";
$tools{sort2temp_pl}      = "/home/wangming/work/bin/get_sRNA/sort2temp.pl";
$tools{soapPE_merge_pl}   = "/home/wangming/work/bin/get_sRNA/soapPE_merge.pl";
$tools{soap2sam_pl}       = "/home/wangming/work/bin/soap2sam.pl";
$tools{perl}              = "/home/wangming/localperl/bin/perl";
$tools{htseq_count}       = "/home/wangming/.local/bin/htseq-count";
$tools{bedtools}          = "/home/wangming/software/bedtools2-2.21.0/bin/bedtools";
$tools{sort2bed_pl}       = "/home/wangming/work/bin/sort2bed.pl";
&Check_file_exists(values %tools);

################################################################################
# Step 1. Prepare reads and reference files
die "File name should match: 1.[a-zA-Z0-9_], 2.End with a letter[a-zA-Z]; eg: H37Rv, B42C\n" unless ($opts{n} =~ /^\w+\D$/);

# print "\# Check 1.Reference.fa/.gff/.fa.fai files\n";
my $Check_ref = $opts{d};
my $Check_gff = $opts{f};
my $Check_fai = $opts{d}.".fai";
#&Check_file_exists($Check_ref, $Check_gff, $Check_fai);   # Check file existence

# 1). Check reference file
die "Ref file Not found:\n$opts{d}" unless -e $opts{d};
die "Ref file Not found:\n$opts{f}" unless -e $opts{f};
my $Ref_genome = abs_path($opts{d});
my $Ref_genome_path = dirname($Ref_genome);
my $Ref_genome_fai = $Ref_genome.".fai";
system "$tools{samtools} faidx $Ref_genome " unless -e $Ref_genome_fai;

# 2). Check Index files
my $Aligner_index;
my $Aligner_index_bowtie_check = $Ref_genome_path."/".$opts{n}.".1.ebwt";
my $Aligner_index_bowtie2_check = $Ref_genome_path."/".$opts{n}.".1.bt2";
my $Aligner_index_soap_check   = $Ref_genome.".index.amb";

if($opts{p} eq "bowtie"){
    system"$tools{Aligner_build} -q -f $Ref_genome $Ref_genome_path\/$opts{n}" unless -e $Aligner_index_bowtie_check;
    $Aligner_index = $Ref_genome_path."/".$opts{n};
# print "bowtie\n";
}elsif($opts{p} eq "bowtie2"){ 
    system"$tools{Aligner_build} -q -f $Ref_genome $Ref_genome_path\/$opts{n}" unless -e $Aligner_index_bowtie2_check;
    $Aligner_index = $Ref_genome_path."/".$opts{n};
# print "bowtie2\n";
}elsif($opts{p} eq "soap"){
    system "$tools{Aligner_build} $Ref_genome" unless -e $Aligner_index_soap_check;
    $Aligner_index = $Ref_genome.".index";
#print "soap\n";
}else{
    die "Input -p should be bowtie/bowtie2/soap";
}

# 3). Search read files
my $reads_abs_path = abs_path($opts{r});
my @reads = glob"$opts{r}\/$opts{n}*.f*";
die "File Not Found:\n$reads_abs_path\/$opts{n}01.fa or $opts{n}03_1.fastq" if(@reads < 1);
my @reads_libs = &Get_reads_libs(@reads);
## Find unique lib names
my %tmp = ( );
foreach (@reads_libs){$tmp{$_}++; }
my @reads_libs_uniq = keys %tmp;

sub Get_reads_libs{
    my @file_libs;
    foreach(@_) {
        die "File Not Fount:\n$_" unless -e $_;
        my $file_name  = basename($_);
        my ($file_lib) = $_ =~ /$opts{n}(\d+)/;
        push @file_libs, $file_lib;
    }
    return @file_libs;
}

################################################################################
# Generate command lines

my $pTime = '
Date=$(date)
echo "$Date  ';

my @Align_cmd_lines = ();
push @Align_cmd_lines, "\#\!\/bin\/bash\n";
push @Align_cmd_lines, "$pTime \t \# Step 1. Filter reads by quality\" \n";

################################################################################
# Step 1. Filter reads by quality

########################## Need improvement for multiple samples################
# Step 2. Create Aligner command lines
#my @libs = ("01", "02", "03", "04"); 
my @libs = sort {$a<=>$b} @reads_libs_uniq;
# 01=18-40nt, SE; 02=40-80nt, SE; 03=80-140nt, PE; 04=>140nt, PE
## Create 02.Alignment directory at current work dir if it does not exist
mkdir $opts{A} unless -d $opts{A};
my $work_dir = `pwd`; chomp($work_dir);
my $Alignment_path = $work_dir."/".$opts{A};

push @Align_cmd_lines, "$pTime \t \# Step 2. Mapping reads to reference \" \n";

if($opts{p} eq "soap") {
###############################################################################
# SOAP commands
    my $soap_cmd_SE = " -v 2 -r 2 -p 2 -M 4 ";
    my $soap_cmd_PE = " -m 2 -r 2 -p 2 -x 500 -s 40 -l 32 ";
    foreach my $lib (@libs){
##############################
# A. Create SOAP command lines
        push @Align_cmd_lines, "\n$pTime \t \#\#\# Mapping: $opts{n}$lib \" ";
        my @Reads2bam_cmd = ();
        my $lib_out_prefix = $Alignment_path."/".$opts{n}.$lib;
        if($lib eq "01" || $lib eq "02"){
            my @reads = <$reads_abs_path\/$opts{n}$lib\.fa*>;  # search SE fa/fastq
            die "Find more $reads_abs_path\/$opts{n}\.$lib\. files " if(@reads > 1); ###
            my $lib_reads = shift(@reads);
#            my $lib_reads = $reads_abs_path."/".$opts{n}.$lib.".fa";  ##
#            my $lib_out_prefix = $Alignment_path."/".$opts{n}.$lib;
            my $R_soap = $lib_out_prefix.".soap";
            my $R_unsoap = $lib_out_prefix.".unsoap";        
            $Reads2bam_cmd[0] = join" ",($tools{Aligner}, $soap_cmd_SE, "-D", $Aligner_index, "-a", $lib_reads, "-o", $R_soap, "-u", $R_unsoap );
            &Check_file_exists($tools{Aligner}, $lib_reads);
##############################
# 1). SE soap -> sam
            my $soapSE_sam = $lib_out_prefix.".sam";
            my $soapSE_bam = $lib_out_prefix.".bam";
            my $soapSE_bam_sorted = $lib_out_prefix.".s.bam";
            $Reads2bam_cmd[1] = "$tools{perl} $tools{soap2sam_pl}  $R_soap  > $soapSE_sam";
            $Reads2bam_cmd[2] = "$tools{samtools} view -bt $Ref_genome\.fai $soapSE_sam -o $soapSE_bam";
            $Reads2bam_cmd[3] = "$tools{samtools} sort $soapSE_bam $lib_out_prefix\.s";
            $Reads2bam_cmd[4] = "$tools{samtools} index $soapSE_bam_sorted";
            push @Align_cmd_lines, @Reads2bam_cmd;
#        }elsif($lib eq "03" || $lib eq "04") {
        }elsif($lib > "02"){
            my $lib_reads_1 = $reads_abs_path."/".$opts{n}.$lib."_1.fastq";
            my $lib_reads_2 = $reads_abs_path."/".$opts{n}.$lib."_2.fastq";
#            my $lib_out_prefix = $Alignment_path."/".$opts{n}.$lib;
            my $R_soap = $lib_out_prefix.".PEsoap";
            my $R_unsoap = $lib_out_prefix.".PEunsoap";
            my $R_Singlesoap = $lib_out_prefix.".Singlesoap";
            $Reads2bam_cmd[0] = join" ",($tools{Aligner}, $soap_cmd_PE, "-D", $Aligner_index, "-a", $lib_reads_1, "-b", $lib_reads_2, "-o", $R_soap, "-2", $R_Singlesoap, "-u", $R_unsoap);
            &Check_file_exists($tools{Aligner}, $lib_reads_1, $lib_reads_2);
##############################
# 2). PE soap -> sam
            my $soapPE_merge = $lib_out_prefix.".soap";
            my $soapPE_merge_fa = $lib_out_prefix.".fa";
            my $soapPE_merge_del = $lib_out_prefix.".del";
            my $soapPE_merge_sam = $lib_out_prefix.".sam";
            my $soapPE_merge_bam = $lib_out_prefix.".bam";
            my $soapPE_merge_bam_sorted = $lib_out_prefix.".s.bam";
            $Reads2bam_cmd[1] = "$tools{perl} $tools{soapPE_merge_pl} -a $Ref_genome -i $R_soap -o $soapPE_merge -c $soapPE_merge_fa -d $soapPE_merge_del";
            $Reads2bam_cmd[2] = "$tools{perl} $tools{soap2sam_pl} $soapPE_merge  > $soapPE_merge_sam";
            $Reads2bam_cmd[3] = "$tools{samtools} view -bt $Ref_genome\.fai $soapPE_merge_sam -o $soapPE_merge_bam";
            $Reads2bam_cmd[4] = "$tools{samtools} sort $soapPE_merge_bam  $lib_out_prefix\.s ";
            $Reads2bam_cmd[5] = "$tools{samtools} index $soapPE_merge_bam_sorted";
            push @Align_cmd_lines, @Reads2bam_cmd;
        }else{
            die "Check the value of @libs, should be 01-04 ";
        }
#        push @Align_cmd_lines,"\n";
    }
} elsif($opts{p} eq "bowtie") {
###############################################################################
# B. Create Bowtie commands
    my $bowtie_cmd_SE = " -v 3 --best --strata -S -m 100 -X 300 --chunkmbs 256 -p 4 ";
    my $bowtie_cmd_PE = " -v 3 --best --strata -S -m 100 -I 60 -X 300 --chunkmbs 256 -p 4 ";
    # SE  bowtie -para index reads
    foreach my $lib (@libs){
        my $lib_out_prefix = $Alignment_path."/".$opts{n}.$lib;
        my $Out_sam = $lib_out_prefix.".sam";
        push @Align_cmd_lines,"\n\n$pTime \t \#\#\# Mapping: $opts{n}$lib\" ";
        if($lib eq "01" || $lib eq "02"){
##############################
# 1). SE => sam
            my @reads = <$reads_abs_path\/$opts{n}$lib\.fa*>;
            die "More than one files found: $reads_abs_path\/$opts{n}\.$lib\.fa* " if(@reads > 1);
            my $lib_reads = shift(@reads);
            my $reads_para = "-f";
            $reads_para = "-q" if($lib_reads =~ /\.fastq$/);
            my $SE_cmd = join " ", ($tools{Aligner}, $bowtie_cmd_SE, $Aligner_index, $reads_para, $lib_reads, $Out_sam);
#            my $lib_reads = $reads_abs_path."/".$opts{n}.$lib.".fa";
#            my $SE_cmd = join" ", ($tools{Aligner}, $bowtie_cmd_SE, $Aligner_index, "-f", $lib_reads, $Out_sam);
            &Check_file_exists($tools{Aligner}, $lib_reads);
            push @Align_cmd_lines, $SE_cmd;
#        }elsif($lib eq "03" || $lib eq "04"){
        }elsif($lib > "02"){
##############################
# 2). PE => sam
            my $lib_reads_1 = $reads_abs_path."/".$opts{n}.$lib."_1.fastq";
            my $lib_reads_2 = $reads_abs_path."/".$opts{n}.$lib."_2.fastq";
            my $PE_cmd = join" ", ($tools{Aligner}, $bowtie_cmd_PE, $Aligner_index, "-1", $lib_reads_1, "-2", $lib_reads_2, $Out_sam);
            &Check_file_exists($tools{Aligner}, $lib_reads_1, $lib_reads_2);
            push @Align_cmd_lines, $PE_cmd;
        }else{
            die "Check the value of @libs, should be 01-04 ";
        }
        my $Out_bam = $lib_out_prefix.".bam";
        my $Out_bam_sorted = $lib_out_prefix.".s.bam";
        my $sam_cmd_1 = "$tools{samtools} view -bt $Ref_genome\.fai $Out_sam -o $Out_bam";
        my $sam_cmd_2 = "$tools{samtools} sort $Out_bam $lib_out_prefix\.s";
        my $sam_cmd_3 = "$tools{samtools} index $Out_bam_sorted";
        push @Align_cmd_lines, ($sam_cmd_1, $sam_cmd_2, $sam_cmd_3);
        push @Align_cmd_lines,"\n";
    }
}elsif($opts{p} eq "bowtie2"){
################################################################################
# C. Create Bowtie2 commands
    my $bowtie2_cmd_SE = "-p 4";
    my $bowtie2_cmd_PE = "-I 50 -X 200 -p 4";
    foreach my $lib (@libs){
        my $lib_out_prefix = $Alignment_path."/".$opts{n}.$lib;
        my $Out_sam = $lib_out_prefix.".sam";
        push @Align_cmd_lines, "\n\n$pTime \t \#\#\# Mapping: $opts{n}$lib\" ";
        if($lib eq "01" || $lib eq "02"){
##############################
# 1). SE => sam
#            my $lib_reads = $reads_abs_path."/".$opts{n}.$lib.".fa";
#            my $SE_cmd = join" ", ($tools{Aligner}, $bowtie2_cmd_SE, "-x", $Aligner_index, "-f", $lib_reads, "-S", $Out_sam);
            my @reads = <$reads_abs_path\/$opts{n}$lib\.fa*>;
            die "More than one files found: $reads_abs_path\/$opts{n}\.$lib\.fa* " if(@reads > 1);
            my $lib_reads = shift(@reads);
            my $reads_para = "-f";
            $reads_para = "-q" if($lib_reads =~ /\.fastq$/);
            my $SE_cmd = join " ", ($tools{Aligner}, $bowtie2_cmd_SE, "-x", $Aligner_index, $reads_para, $lib_reads, "-S", $Out_sam);
            &Check_file_exists($tools{Aligner}, $lib_reads);
            push @Align_cmd_lines, $SE_cmd;
#        }elsif($lib eq "03" || $lib eq "04"){
        }elsif($lib > "02"){
##############################
# 2). PE => sam            
            my $lib_reads_1 = $reads_abs_path."/".$opts{n}.$lib."_1.fastq";
            my $lib_reads_2 = $reads_abs_path."/".$opts{n}.$lib."_2.fastq";
            my $PE_cmd = join" ", ($tools{Aligner}, $bowtie2_cmd_PE, "-x", $Aligner_index, "-1", $lib_reads_1, "-2", $lib_reads_2, "-S", $Out_sam);
            &Check_file_exists($tools{Aligner}, $lib_reads_1, $lib_reads_2);
            push @Align_cmd_lines, $PE_cmd;
        }else{
            die "Check the value of @libs, should be 01-04 ";
        }
        my $Out_bam = $lib_out_prefix.".bam";
        my $Out_bam_sorted = $lib_out_prefix.".s.bam";
        my $sam_cmd_1 = "$tools{samtools} view -bt $Ref_genome\.fai $Out_sam -o $Out_bam";
        my $sam_cmd_2 = "$tools{samtools} sort $Out_bam $lib_out_prefix\.s";
        my $sam_cmd_3 = "$tools{samtools} index $Out_bam_sorted";
        push @Align_cmd_lines, ($sam_cmd_1, $sam_cmd_2, $sam_cmd_3);
        push @Align_cmd_lines,"\n";
    }
} else {
    die "Check \$opts{p}, should be soap or bowtie";
}
#print join"\n",@Align_cmd_lines,"\n";

sub Check_file_exists{
    foreach(@_){
        die "File Not Found:\n$_" unless -e $_;
    }
}

################################################################################
# Step 3. BAM => coverage => sRNA
# Create Coverage dir
my $Coverage_dir = $opts{C};
mkdir $Coverage_dir unless -d $Coverage_dir;
my $Coverage_path = abs_path($Coverage_dir);

# BAM to Coverage
push @Align_cmd_lines, "\n$pTime \t \# Step 3. Transfer BAM to Coverage file \" ";
foreach my $lib (@libs) {
    my $bam_in = $work_dir."/".$opts{A}."/".$opts{n}.$lib.".s.bam";
#    &Check_file_exists($bam_in);
    my $bam_in_abs_path = abs_path($bam_in);
    my $bam_name = basename($bam_in);
    my $Coverage_prefix = $opts{n}.$lib;
    my @bam2cov_cmd = ();
    $bam2cov_cmd[0] = "\n\# Sample: $opts{n}$lib ";
    $bam2cov_cmd[1] = "$tools{genomeCoverageBed} -ibam $bam_in_abs_path -d -g $opts{n}\.fai -strand - > $Coverage_path\/$Coverage_prefix\.coverage\.n ";
    $bam2cov_cmd[2] = "$tools{genomeCoverageBed} -ibam $bam_in_abs_path -d -g $opts{n}\.fai -strand + > $Coverage_path\/$Coverage_prefix\.coverage\.p ";
    $bam2cov_cmd[3] = "$tools{perl} $tools{get_sRNA_pl} -i $Coverage_path\/$Coverage_prefix\.coverage\.n -c 0.5 -s - -o $Coverage_path\/$Coverage_prefix\_tag.tmp\.n ";
    $bam2cov_cmd[4] = "$tools{perl} $tools{get_sRNA_pl} -i $Coverage_path\/$Coverage_prefix\.coverage\.p -c 0.5 -s + -o $Coverage_path\/$Coverage_prefix\_tag.tmp\.p ";
    $bam2cov_cmd[5] = "cat $Coverage_path\/$Coverage_prefix\_tag.tmp\.n $Coverage_path\/$Coverage_prefix\_tag.tmp\.p > $Coverage_path\/$Coverage_prefix\_tag.tmp ";
    $bam2cov_cmd[6] = "$tools{perl} $tools{getPosition_pl} -i $Coverage_path\/$Coverage_prefix\_tag.tmp -f $opts{f} -o $Coverage_path\/$Coverage_prefix\_tag\.txt ";
    $bam2cov_cmd[7] = "$tools{perl} $tools{sort2candi_pl}  $Coverage_path\/$Coverage_prefix\_tag\.txt ";
    push @Align_cmd_lines, @bam2cov_cmd;
}

################################################################################
# Step 4. Merge sRNAs from all libraries
mkdir $opts{M} unless -d $opts{M}; # 04.Merge_tags
my $merge_path = abs_path($opts{M});

push @Align_cmd_lines, "\n$pTime \t \# Step 4. Merge sRNA tags from all libraries\" \n";
my @tag_files;
foreach my $lib (@libs){
    my $t_file = "$Coverage_path\/$opts{n}$lib\_tag\.txt";
    push @tag_files, $t_file;
}
my $parameter_file = join" ",(@tag_files);
my $merge_cmd = "$tools{perl} $tools{merge4all_pl}  $parameter_file  > $merge_path\/$opts{n}\_merge.tmp ";
my $merge2temp = "$tools{perl} $tools{sort2temp_pl} $merge_path\/$opts{n}\_merge.tmp";
my $merge2temp2pos = "$tools{perl} $tools{getPosition_pl}  -i $merge_path\/$opts{n}\_merge.tmp.temp -f $opts{f} -o $merge_path\/$opts{n}\_merge.txt";
my $merge2sort = "$tools{perl} $tools{sort2candi_pl}  $merge_path\/$opts{n}\_merge.txt";
push @Align_cmd_lines, ($merge_cmd, $merge2temp, $merge2temp2pos, $merge2sort);

push @Align_cmd_lines, "\n$pTime \t \# Finsh. Write sRNAs to 04.Merge/$opts{n}\_merge.txt \" ";

#print join"\n",@Align_cmd_lines,"\n";

################################################################################
# Step 5. Count reads on each sRNAs
# 1. htset-count -f bam  --stranded=yes -q Alignment.s.bam in.gff > out.exp
# 2. bedtools multicov -s -bams align1.bam align2.bam ... -bed in.gff >out.exp
push @Align_cmd_lines, "$pTime \t \# Step 5. Count reads on each sRNAs\" ";

mkdir $opts{E} unless -d $opts{E};
my $count_path = abs_path($opts{E});

my @Step5cmd = ();
my $sRNA_file = abs_path("$merge_path\/$opts{n}_merge\.txt");
my $sRNA_gff  = basename($sRNA_file);
$sRNA_gff =~ s/\.txt/.gff/;
$sRNA_gff = $count_path."/".$sRNA_gff;
my $trans_sort2gff = "$tools{perl} $tools{sort2bed_pl} -t sort2gff -i $sRNA_file -o $sRNA_gff  -g $Ref_genome ";
my $extra_1 = "head $sRNA_gff > tmp.gff";
my $extra_2 = "mv -f tmp.gff $sRNA_gff";
push @Align_cmd_lines, $trans_sort2gff;
push @Align_cmd_lines,"\n";
#push @Step5cmd, ($trans_sort2gff, $extra_1, $extra_2, "\n");

#&Check_file_exists($sRNA_file); ## Update Nov-17

foreach my $lib (@libs){
    my $sorted_bam = "$Alignment_path\/$opts{n}$lib\.s\.bam";
#    warn "BAM file Not Found:\n$sorted_bam" unless -e $sorted_bam; ## Update Nov-17
    my $lib_cmd = &CountReads($sRNA_file, $sorted_bam, $lib);
    push @Align_cmd_lines, "$pTime \t \# Count lib: $lib \" ";
    push @Align_cmd_lines, $lib_cmd;
#    push @Step5cmd, ("$pTime  \# Run lib: $lib \" ", $lib_cmd);
}

#print join"\n",@Step5cmd;
print join"\n", @Align_cmd_lines, "\n";

sub CountReads{
    my ($in_txt, $in_bam, $lib_id) = @_;
    $in_txt = abs_path($in_txt);
    $in_bam = abs_path($in_bam);
# Prepare output file names    
    my $in_txt_exp = basename($in_txt);
    $in_txt_exp =~ s/txt$/lib$lib_id\.txt/;
    $in_txt_exp = $count_path."/".$in_txt_exp;
# Create commands
    my $multicov_cmd = "$tools{bedtools} multicov -s -bams $in_bam -bed $sRNA_gff > $in_txt_exp ";
    return join"\n",($multicov_cmd, "\n");
#    my $htseq_cmd   = "$tools{htseq_count}  -f bam -q -o out.sam  $in_bam  $sRNA_gff  > $in_txt_exp";    
#    return join"\n",($htseq_cmd,"\n");  
}

################################################################################
# Create PBS template
my $PBS_Header = "#PBS -S /bin/bash
#PBS -q Debug
#PBS -o $work_dir/pbs.out
#PBS -e $work_dir/pbs.err
#
cd \$PBS_O_WORKDIR
#nohup  sh YOUR.sh > YOUR.sh.log &

# Write your Command here ... : Check the queue: Debug/General
# sh  \$PBS_O_WORKDIR/work.sh       > \$PBS_O_WORKDIR/work.sh.log
#
";

open H,"> wm_job.pbs" or die;
print H $PBS_Header;

