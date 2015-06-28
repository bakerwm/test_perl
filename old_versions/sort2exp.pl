### !---------------------------------------------
### delete, replaced by: chk_feature2count.pl

#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use File::Basename qw/basename/;

################################################
# This script is designed to count the number 
# of reads mapping to sRNA/genes, using BAM files
#
# Tools: BEDTools/coverageBed
# Date: 2014-01-14
################################################

my %opts = ();
getopts("a:i:o:",\%opts);

die "Usage: perl $0 -a in.bam -i sRNA.txt -o out.count.txt" unless($opts{a} && $opts{i});

# Tools
my $coverageBed = '/home/wangming/Documents/bedtools-2.17.0/bin/coverageBed';
my $samtools = '/home/wangming/Documents/samtools-0.1.19/samtools';

# Input files
my $InName = basename($opts{i});
my $OutBED = my $OutCOUNT = $InName; 
$OutBED =~ s/\.txt/.bed/;
$OutCOUNT =~ s/\.txt/.count.txt/;
$opts{o} = $OutCOUNT unless($opts{o});

# Get name of chromosome/reference
my $chr = "chr";
my @heads = `$samtools view -H $opts{a} `;  # Run the samtools view
foreach(@heads){
    next unless(/^\@SQ/);
    ($chr) = /SN\:(.*)\sLN\:/;  # @SQ     SN:NC_000913.2  LN:4639675
}

# Trans Infile to BED format;
open F,"$opts{i}" or die "Cannot open $opts{i}, $!";
my @RAs = <F>;
open O,"> $OutBED" or die "Cannot open $OutBED, $!";
foreach(@RAs){
    chomp;
    next if((split/\t/) < 6);
    my $out = join"\t",($chr,(split/\t/)[3,4,0,2,5]);
    next if((split/\t/,$out)<6);
    print O join"\t",$out,"\n";
}
close F;
close O;

# Run coverageBed
my @count = `$coverageBed -abam $opts{a} -b $OutBED -s`; # Run coverageBed, count reads on sRNAs.

# Get the exp result
my %COUNT = ();
foreach(@count){
    my($id,$exp) = (split/\t/)[3,6];
    $COUNT{$id} = $exp;
}

# Write the exp into InTxt file
open OUT,"> $opts{o}" or die "Cannot open file $opts{o}, $!";
foreach(@RAs){
    chomp;
    my @tabs = split/\t/;
    next if(@tabs < 6 || /^\s/);
    my $p = $tabs[1];
#    die "Check the Count of $tabs[0] " unless($COUNT{$tabs[0]}); ## Check the expression
    $tabs[1] = $COUNT{$tabs[0]};
    print OUT join "\t",(@tabs,$p),"\n";
}
close OUT;

#unlink ("$OutBED");

