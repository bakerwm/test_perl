#!/usr/local/bin/perl -w

#################################################
# Convert BAM file to bedgraph 
#
#################################################

use strict;
use warnings;
use File::Basename qw/basename dirname/;

sub help{<<EOF
Usage: perl  Bam2bedgraph.pl  ref.fa ref.list  ref.bed  sorted.bam/
    ref.fna:    begin as >gi|49175990|ref|NC_000913.2| or >NC_000913.2
    ref.list:   as genome list "NC_000913.2  4639675"
    ref.bed:    as genome length "NC_000913.2 0 4639675 Ecoli 1 +"
    bam_dir:    as "Alignment/"
    
### Results output to the same directory of bam files.
        
Example:
    perl Bam2bedgraph.pl NC_000913.fna ref.list  ref.bed  Alignment/

EOF

}

my $ref_fa    = shift || die &help;
my $ref_list  = shift || die &help;
my $ref_bed   = shift || die &help;  
my $bam_dir   = shift || die &help;

## Tools
my $genomeCoverageBed = '/home/wangming/Documents/bedtools-2.17.0/bin/genomeCoverageBed';
my $coverageBed       = '/home/wangming/Documents/bedtools-2.17.0/bin/coverageBed';

## Find bam files
my @SortedBam = <$bam_dir*sorted.bam>;
my $num = @SortedBam;
print "Find $num sorted Bam files:\n";

## Find gene id
open FA, "$ref_fa" or die "Cannot open $ref_fa $! \n";
my $ref_head = <FA>;
my ($ref_gi,$ref_gi_id);
if($ref_head =~ /^\>gi\|\d+\|ref/){
    ($ref_gi) =  $ref_head =~ /\>(gi\|\d+\|ref\|.*\|)\s/;    # >gi|49175990|ref|NC_000913.2| Escherichia ...
    $ref_gi_id = (split/\|/,$ref_gi)[3];
}
close FA;

### Check the report file
my $report = dirname($SortedBam[0])."/Genome_map.report";
unlink $report if -e $report;

## output bedgraph files & coverage report;
foreach my $i (@SortedBam){
print "$i\n";
    my $fname = $i;
    $fname =~ s/\.sorted.bam//;
    my $bedgraphN = $fname.".Neg.bedgraph";
    my $bedgraphP = $fname.".Pos.bedgraph";

## Create BedGraph files
    if($ref_head =~ /^\>gi\|\d+\|ref/){
        system `$genomeCoverageBed -ibam $i -bga -g $ref_list -strand - | sed -e 's/$ref_gi/$ref_gi_id/g' > $bedgraphN`;
        system `$genomeCoverageBed -ibam $i -bga -g $ref_list -strand + | sed -e 's/$ref_gi/$ref_gi_id/g' > $bedgraphP`;
    }

    system `$genomeCoverageBed -ibam $i -bga -g $ref_list -strand - > $bedgraphN`;
    system `$genomeCoverageBed -ibam $i -bga -g $ref_list -strand + > $bedgraphP`;

### Create genome coverage report
    system `$coverageBed  -abam $i -b $ref_bed -s >> $report `; 
}
    
