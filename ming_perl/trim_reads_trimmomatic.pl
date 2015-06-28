#!/home/wangming/localperl/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use Cwd 'abs_path';
use File::Basename qw/basename dirname/;

sub usage{
    print << "EOF"; 
Trim the RNA-Seq reads and Align to reference genome using bowtie. And statistic the
alignment result.

Note:
    -i :    Input the fastq file of SE or -1 of PE
    -b :    The -2 reads files of PE
    -t :    SE/PE
    -o :    The output of commands

example:
1. For SE reads
perl TrimReads.pl -i input.fastq -t SE -o command.sh

2. For PE reads
perl TrimReads.pl -i input_1.fastq -b input_2.fastq -t PE -o command.sh

EOF
exit(1);
}

my %opts = ();
getopt("i:b:t:o:", \%opts);
&usage unless (defined $opts{i} && defined $opts{t});

################################################################################
# Tools
my $trim_jar_path = "/home/wangming/software/Trimmomatic-0.32";
my $trim_adapters = $trim_jar_path."\/adapters\/TruSeq2";
my $bowtie    = '/home/wangming/software/bowtie-1.1.1/bowtie';
my $H37Rv_build = '/home/wangming/work/database/H37Rv/H37Rv';
my $samstat   = '/home/wangming/software/samstat-1.5/bin/samstat';
my $java      = '/usr/bin/java';

my $bowtie_para_SE = "-v 3 --best --strata -S -m 100 -X 300 --chunkmbs 256 -p 4";
my $bowtie_para_PE = "-v 3 --best --strata -S -m 100 -I 60 -X 300 --chunkmbs 256 -p 4";

my $trim_para_SE = "ILLUMINACLIP\:$trim_adapters\-SE.fa\:2\:30\:10 LEADING\:3 TRAILING\:3 SLIDINGWINDOW\:4:15";
my $trim_para_PE = "ILLUMINACLIP\:$trim_adapters\-PE.fa\:2\:30\:10 LEADING\:3 TRAILING\:3 SLIDINGWINDOW\:4\:15";


################################################################################
# 1. Trim reads
# 2. bowtie, sam2bam
# 3. samstat *.sam *.fastq
# 4. samtools flagstat

my $work_dir = `pwd`; chomp($work_dir);
mkdir "$work_dir\/reads_trim" unless -d "$work_dir\/reads_trim";
my $trim_path = $work_dir."\/reads_trim";

my @cmd_out = ();

push @cmd_out,("\#/bin/bash\n\n");

# For SE reads
if($opts{t} eq "SE") {
    die "File Not Found: $opts{i}" unless -e $opts{i};
    my $in_reads = abs_path($opts{i});
    my $trim_reads = basename($opts{i});
    die "Need input Fastq file." unless ($opts{i} =~ /\.fastq$/);
    $trim_reads =~ s/\.fastq$/.trim.fastq/;
    $trim_reads = $work_dir."\/reads_trim\/$trim_reads";
    my $cmd_trim = "$java -jar $trim_jar_path\/trimmomatic-0.32.jar SE $in_reads $trim_reads $trim_para_SE MINLEN\:15 ";

    my ($bowtie_out) = basename($opts{i})=~ /^(.*)\.+.*\.fastq$/;
    my $bowtie_out_raw = $trim_path."\/".$bowtie_out.".raw.sam";
    my $bowtie_out_trim = $trim_path."/".$bowtie_out.".trim.sam";
    my $cmd_bowtie1 = "$bowtie $bowtie_para_SE $H37Rv_build -q $in_reads $bowtie_out_raw ";
    my $cmd_bowtie2 = "$bowtie $bowtie_para_SE $H37Rv_build -q $trim_reads $bowtie_out_trim ";
    my $cmd_stat1 = "$samstat $bowtie_out_raw";
    my $cmd_stat2 = "$samstat $bowtie_out_trim";
    push @cmd_out, ($cmd_trim, $cmd_bowtie1, $cmd_bowtie2, $cmd_stat1, $cmd_stat2);
# For PE reads
} elsif($opts{t} eq "PE" ) {
    die "Need -b read2.fastq for PE." unless(defined $opts{b});
    die "File Not Found: $opts{i}" unless -e $opts{i};
    die "File Not Found: $opts{b}" unless -e $opts{b};
    my $in_a_reads = abs_path($opts{i});
    my $in_b_reads = abs_path($opts{b});
    my $trim_a_name = basename($opts{i}); $trim_a_name =~ s/\.fastq$//; 
    my $trim_b_name = basename($opts{b}); $trim_b_name =~ s/\.fastq$//;
    my $trim_a_paired = $trim_path."/".$trim_a_name.".paired.fastq";
    my $trim_a_unpaired = $trim_path."/".$trim_a_name.".unpaired.fastq";
    my $trim_b_paired = $trim_path."/".$trim_b_name.".paired.fastq";
    my $trim_b_unpaired = $trim_path."/".$trim_b_name.".unpaired.fastq";
    my $cmd_trim = "$java -jar $trim_jar_path\/trimmomatic-0.32.jar PE $in_a_reads $in_b_reads $trim_a_paired $trim_a_unpaired $trim_b_paired $trim_b_unpaired $trim_para_PE MINLEN\:40";
 
    my ($bowtie_out) = basename($opts{i}) =~ /(^\w+)\_\d.*\.fastq$/;
    my $bowtie_out_raw = $trim_path."\/".$bowtie_out.".raw.sam";
    my $bowtie_out_trim = $trim_path."\/".$bowtie_out.".trim.sam";
    my $cmd_bowtie1 = "$bowtie $bowtie_para_PE $H37Rv_build -1 $in_a_reads -2 $in_b_reads $bowtie_out_raw";
    my $cmd_bowtie2 = "$bowtie $bowtie_para_PE $H37Rv_build -1 $trim_a_paired -2 $trim_b_paired $bowtie_out_trim";
    my $cmd_stat1 = "$samstat $bowtie_out_raw";
    my $cmd_stat2 = "$samstat $bowtie_out_trim";
    push @cmd_out, ($cmd_trim, $cmd_bowtie1, $cmd_bowtie2, $cmd_stat1, $cmd_stat2);    
}

print join"\n", @cmd_out,"\n";



my $pbs_header = "#PBS -S /bin/bash
#PBS -q Debug
#PBS -o $work_dir\/pbs.out
#PBS -e $work_dir\/pbs.err

# Write your Command here ... : Check the queue: Debug/General
# sh  $work_dir\/bowtie_align.sh       > $work_dir\/bowtie_align.sh.log
 
";

open OUT,"> test.pbs" ;
print OUT $pbs_header;
close OUT;
