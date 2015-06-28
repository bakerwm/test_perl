#!/usr/bin/perl  -w
use strict;

my $usage = "perl  get_sRNA_sh.pl  H37Rv
Description:
This script is designed to get sRNAs.
Create two directories:
1. GFF file.    in dir: database/
2. sRNAs in each libraries.  in dir: raw_data/
";

my $strain = shift or die $usage;

##################
# Perl scripts
##################
my $merge4all  = '/home/wangming/work/bin/get_sRNA/merge4all.pl';
my $sort2temp  = '/home/wangming/work/bin/get_sRNA/sort2temp.pl';
my $sort2candi = '/home/wangming/work/bin/get_sRNA/sort2candi_v1.pl';
my $getPosition = '/home/wangming/work/bin/get_sRNA/getPosition.pl';

##################
# Search rpkm files
##################
my @sRNA_rpkm = glob"raw_data/*.txt";
my $infile    = join"\,", @sRNA_rpkm;
my $number    = @sRNA_rpkm;
my $GFF       = glob"database/*gff";

system"perl  $merge4all  -n $number  -i $infile  -o  temp.txt";
system"perl  $merge4all  -n 1        -i temp.txt -o  temp_v1.txt";
system"perl  $merge4all  -n 1        -i temp_v1.txt -o temp_v2.txt";
system"perl  $merge4all  -n 1        -i temp_v2.txt -o temp_v3.txt";
system"perl  $merge4all  -n 1        -i temp_v3.txt -o temp_v4.txt";
system"perl  $merge4all  -n 1        -i temp_v4.txt -o temp_v5.txt";
system"mv temp_v5.txt  $strain\_merge\.txt";
unlink glob "temp_v*";
system"perl $sort2temp  $strain\_merge\.txt";
system"perl $getPosition $GFF $strain\_merge.txt.temp  $strain  >$strain\_merge_lncRNA.txt";
system"perl $sort2candi  $strain\_merge_lncRNA.txt";

if(-d "$strain\_sRNA"){
    system"rm -r  $strain\_sRNA";
}
#else{
#    mkdir ("$strain\_sRNA", 0755);
#}
#system"mv  $strain\_merge*  temp*   $strain\_sRNA/ ";
