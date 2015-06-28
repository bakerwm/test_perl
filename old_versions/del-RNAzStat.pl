### !------------------------------------
### replaced by: chk_seq2rnaz.pl

#!/home/wangming/localperl/bin/perl  -w
use strict;
#use warnings;  ## Confict with 

use lib qw(/home/wangming/localperl/share/perl5
           /home/wangming/localperl/lib/5.16.1
           /home/wangming/localperl/lib/site_perl/5.16.1
       );

use Getopt::Std;
#use FindBin;
#use lib "$FindBin::Bin/../lib";
#use File::Temp qw(tempdir);
use File::Basename qw(dirname basename);
use File::Spec::Functions qw(catdir catfile);
use File::Which qw(which where);
use Cwd qw(abs_path cwd);

use Data::Dumper;

sub usage{
    print << "EOF";
Usage: perl RNAzStat.pl -i input.maf -n known.bed

  -i : input MAF multi sequences alignments results, (ClustalW or Clustal Omega, Multiz)
  -n : The file of known ncRNAs in BED format.


EOF
exit;
}

### pipeline ###
# 1. rnazWindow.pl --min-seqs=4  input.maf > windows.maf
# 2. RNAz --both-strands --no-shuffle --cutoff=0.5 windows.maf > rnaz.out
# 3. rnazCluster.pl rnaz.out > results.dat
# 
# Assess the *.dat data
# rnazIndex.pl --bed results.dat | rnazBEDsort.pl | rnazBEDstats.pl > out.stat
#
# Filter
# rnazFilter.pl "P>0.9" results.dat | rnazSort.pl combPerPair > filt.dat
#
# Annotation
# 1. using BED file
# rnazAnnotate.pl --bed known.bed results.dat > annotated.dat
# 2. using Rfam database
# rnazBlast.pl --database rfam --seq-dir=seq --blast-dir=rfam results.dat > annotated.dat
#
# Visualize
# rnazCluster.pl --html rnaz.out > results.dat # create results/ dir,
# rnazIndex.pl --html results.dat > results/results.html # summary html file.
#
# Estimate False Positives
# 1. Get random samples (shuffle input.maf)
# rnazRandomizeAln.pl input.maf > random-input.maf
# 2. Run RNAz and statistic procedure
# rnazWindow.pl  | RNAz | rnazCluster.pl | rnazFilter.pl | rnazAnnotate.pl | rnazBEDstats.pl | 
#
### Output ###
# 1. annotated.dat, bed, gff file
# 2. stat files, (P>0.5, P>0.8)



my %opts = ();
getopt("i:n:o:",\%opts);
&usage unless($opts{i});

my %Tools = ();
######
# Tools RNAz
my $RNAz_path = which('RNAz');
my $RNAz_perl_dir = catdir(dirname($RNAz_path), '../perl');
my @RNAz_perls = glob("$RNAz_perl_dir\/*.pl");
foreach my $p (@RNAz_perls){
    my $perl_name = basename($p);
    $Tools{$perl_name} = abs_path($p);
}

my $original_dir = cwd;



#######################################
# Create temp dir
mkdir 'temp' unless -d 'temp';
unlink './temp/*';

#######################################
# Check input maf
my $input_maf = abs_path($opts{i});
open F, $input_maf or die "Cannot open $opts{i}, $!";
my $check_maf = <F>;
die "Input file should be MAF format: ##maf version" unless($check_maf =~ /^\#\#maf/);
close F;

# 1. processing the raw alignmnts (MAF)
my $in_ID = basename($opts{i}); $in_ID =~ s/\.maf$//i;
my $out_dir = catdir(cwd, "temp");
my $window_maf = catfile($out_dir, "$in_ID\_windows.maf");
my $rnaz_out = catfile($out_dir, "$in_ID\_rnaz.out");
my $rnaz_dat = catfile($out_dir, "$in_ID.dat");
my $h1 = "\#\# Processing $in_ID\.maf with rnazWindow.pl ";
my $p1 = "rnazWindow.pl --min-seqs=4 $opts{i} > $window_maf";
my $p2 = "RNAz --both-strands --no-shuffle --cutoff=0.5  $window_maf > $rnaz_out";
my $p3 = "rnazCluster.pl $rnaz_out > $rnaz_dat";

# 2. Filtering results
my @p4 = &RNAz_dat("$rnaz_dat");

# 3. Annotate by BED file
my @p5 = &RNAz_annotate("$rnaz_dat");

# 4. Annotate by Blast rfam



## Create PBS job
my $job_header = "\#PBS -S /bin/bash
\#PBS -N RNAz_stat.job
\#PBS -o $original_dir\/RNAz_stat.job.out
\#PBS -e $original_dir\/RNAz_stat.job.err
\#PBS -q General
\#
";

print join"\n", ($job_header,$h1,$p1,$p2,$p3,@p4,@p5), "\n";


sub RNAz_dat{
    my $dat = shift(@_);
    my $out_dir = catdir(cwd, "temp");
    my $dat_name = basename($dat); $dat_name =~ s/\.dat$//g;
    my $stat_raw = catfile($out_dir, "$dat_name\_raw.stat");
    my $stat_P9  = catfile($out_dir, "$dat_name\_P0.9.stat");
    my $stat_P5  = catfile($out_dir, "$dat_name\_P0.5.stat");
    my $stat1 = "rnazIndex.pl --bed $dat \| rnazBEDsort.pl \| rnazBEDstats.pl > $stat_raw";
    my $stat2 = "rnazFilter.pl \"P>=0.5\" $dat \| rnazIndex.pl --bed \| rnazBEDsort.pl \| rnazBEDstats.pl > $stat_P5";
    my $stat3 = "rnazFilter.pl \"P>=0.9\" $dat \| rnazIndex.pl --bed \| rnazBEDsort.pl \| rnazBEDstats.pl > $stat_P9 ";
    return ($stat1, $stat2, $stat3);
}

sub RNAz_annotate{
    my $dat = shift(@_);
#    die "File not exists: $dat \n" unless -e $dat;
    die "Known seq file not found: $opts{n} \n" unless -e $opts{n};
    my $out_dir = catdir(cwd, "temp");
    my $dat_name = basename($dat); $dat_name =~ s/\.dat$//g;
    my $anno_raw = catfile($out_dir, "$dat_name\_raw.anno");
    my $anno_P9 = catfile($out_dir, "$dat_name\_P0.9.anno");
    my $anno_P5 = catfile($out_dir, "$dat_name\_P0.5.anno");
    my $annotate1 = "rnazAnnotate.pl --bed $opts{n} $dat > $anno_raw";
    my $annotate2 = "rnazFilter.pl \"P>0.5\" $dat \| rnazAnnotate.pl --bed $opts{n} > $anno_P5";
    my $annotate3 = "rnazFilter.pl \"P>0.9\" $dat \| rnazAnnotate.pl --bed $opts{n} > $anno_P9";
    return ($annotate1, $annotate2, $annotate3);
}

