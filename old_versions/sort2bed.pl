#!/usr/bin/perl -w
use strict;
use warnings;

use File::Spec::Functions qw(catfile catdir);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path cwd);
use Getopt::Long;
use Pod::Usage;

my %opts = ();
my $type     = '';
my $inFile   = '';
my $outFile   = '';
my $feature  = '';
my $genome   = '';
my $help     = '';
my $man      = '';

GetOptions(
        't|type=s' => \$type,
        'i|input=s' => \$inFile,
        'o|output=s' => \$outFile,
        'f|feature=s' => \$feature,
        'g|genome=s' => \$genome,
        'h|help' => \$help,
        'man' => \$man
        ) or pod2usage(2);

pod2usage(-verbose => 1) if($help);
pod2usage(-verbose => 2) if($man);
pod2usage(-message => "Need choose input file: [-i|--input]") if(! $inFile);

## Choose type of transformer
if(! $type){
    pod2usage(-message => "Need choose transform type: GFF|PTT|BED|Sort, eg: gff2bed");
}else{
    if($type =~ /^(gff|ptt|bed|sort)2(gff|ptt|bed|sort|fa)$/s){
        my ($t1, $t2) = $type =~ /^(\w+)2(\w+)$/;
        die "[-t $type] Choose another type" if ($t1 eq $t2);
    }else{
        pod2usage(-message => "Select [GFF,PTT,BED,Sort], as: gff2bed");
    }
}

## Check input file should be tab-spearated
&CheckFormat($inFile);

##
my ($ta, $tb) = $type =~ /^(\w+)2(\w+)$/;

my %inData = ();
my %outData = ();

# Original format
if($ta =~ /gff/s){
    pod2usage(-message => "Need input [-f|feature] ") if(! $feature);
    %inData = &ReadGFF($inFile, $feature);
}elsif($ta =~ /ptt/s){
    pod2usage(-message => "Need input [-g|genome] ") if(! $genome);
    %inData = &ReadPTT($inFile, $genome);
}elsif($ta =~ /bed/s){
    %inData = &ReadBED($inFile);
}elsif($ta =~ /sort/s){
    %inData = &ReadSort($inFile);
}else{
    
}

$outFile = 'out.txt' if(! $outFile);

open OUT, "> $outFile" or die "$!";
# Target format
if($tb =~ /gff/s){
    my @out = ();
    foreach my $n (sort keys %inData){
       my ($chr, $length, $start, $end, $strand) = split /\t/, $inData{$n};
       my $feature = 'gene';
       my $source = 'RefSeq';
       my $description = "Name=$n;gene=$n;locus_tag=$n";
       my $line = join "\t", ($chr, $source, $feature, $start, $end, '.', $strand, '.', $description);
       push @out, $line;
    }
    my $print_out = join "\n", @out;
    print OUT $print_out, "\n";
}elsif($tb =~ /ptt/s){
    my @out = ();
    foreach my $n (sort keys %inData){
        my ($chr, $length, $start, $end, $strand) = split /\t/, $inData{$n};
        my $colA = $start.'..'.$end;
        my $line = join "\t", ($colA, $strand, '-', '-', $n, $n, '-', '-', $chr);
        push @out, $line;
    }
    my $print_out = join "\n", @out;
    print OUT $print_out, "\n";
}elsif($tb =~ /bed/s){
    my @out = ();
    foreach my $n (sort keys %inData){
        my ($chr, $length, $start, $end, $strand) = split /\t/, $inData{$n};
## 1-st base numbered 0
        $start -= 1;
        $end -= 1;
        my $line = join "\t", ($chr, $start, $end, $n, $length, $strand );
        push @out, $line;
    }
    my $print_out = join "\n", @out;
    print OUT $print_out, "\n";
}elsif($tb =~ /sort/s){
    my @out = ();
    foreach my $n (sort keys %inData){
        my ($chr, $length, $start, $end, $strand) = split /\t/, $inData{$n};
        push @out, join "\t", ($n, $chr, $length, $start, $end, $strand );
    }
    my $print_out = join "\n", @out;
    print OUT $print_out, "\n";
}elsif($tb =~ /fa/s){
    pod2usage(-message => "Need reference file [-g|genome]") if(! $genome);
    my ($id, $ref) = &ReadFa($genome);
    my @out = ();
    foreach my $n (sort keys %inData){
        my ($chr, $length, $start, $end, $strand) = split /\t/, $inData{$n};
        my $seq = &Txt2Seq($ref, $start, $end, $strand);
        push @out, (">$n", $seq);
    }
    my $print_out = join "\n", @out;
    print OUT $print_out, "\n";
}else{

}

####################
# subroutines
####################

## Input should be tab-separated file. [gff/ptt/BED/sort]
sub CheckFormat{
    my $f = shift(@_);
    open F, "< $f" or die "$!";
    while (my $i = <F>){
        if($i =~ /\t/){
            my @tabs = split /\t/, $i;
            die "Input should be tab-separated: $! " if @tabs < 6;
        }
#        die "Input file is not tab-separated: [$f]";
    }
    close F;
}

# Read GFF
sub ReadGFF{
    my ($fh_gff, $feature) = @_;
    my %out = ();
    my $count = 0;
    open F, "<$fh_gff" or die "$!";
    while(my $i = <F>){
        chomp ($i);
        next unless ($i =~ /\t$feature\t/);
        my ($chrom, $start, $end, $strand) = (split /\t/, $i)[0,3,4,6];
        my $length = $end - $start + 1;
        my $name = '';
# feature = gene
        if($feature eq 'gene') {
            if($i=~/locus_tag/){
                ($name) = $i =~ /locus_tag=(\w+)/;
            }elsif(/Name=\w+/){
                ($name) = $i =~ /Name=(\w+)/;
            }else{
                $name = '-';
            }
# feature = CDS|rRNA|tRNA|ncRNA|
        }elsif($feature =~ /CDS|rRNA|tRNA|ncRNA|exon/){
            ($name) = $i =~ /gene=(\w+)/;
# feature not found
        }else{
            pod2usage(-message => "[$feature] may not supported");
        }
        $out{$name} = join "\t", ($chrom, $length, $start, $end, $strand);
        $count ++;
    }
    close F;
    pod2usage(-message => "[$feature] not found in GFF file") if($count == 0);
    return %out;
}

# Read BED
sub ReadBED{
    my $fh_bed = shift(@_);
    my %out = ();
    my $count = 0;
    open F, "<$fh_bed" or die "$!";
    while(my $i = <F>){
        chomp($i);
        next if($i =~ /^\s*$/); # Skip blank lines
        my ($chrom, $start, $end, $name, $strand) = (split /\t/, $i)[0,1,2,3,5];
        $start += 1; ## BED first base = 0;
        $end += 1; ## BED first base = 0;
        my $length = $end - $start + 1;
        $out{$name} = join "\t", ($chrom, $length, $start, $end, $strand);
        $count ++;
    }
    close F;
    pod2usage(-message => "input not like a BED file") if($count == 0);
    return %out;
}

# Read PTT
sub ReadPTT{
    my ($fh_ptt, $g) = @_;
    my %out = ();
    my $count = 0;
    my $chr = '';
    open F, $g or die "Cannot open $g, $!";
    my $line = <F>;
    ($chr) = (split /\|/, $line)[3];
    close F;
    open F2, "<$fh_ptt" or die "$!";
    while(my $i = <F2>){
        chomp ($i);
        next unless($i =~ /^\d+\.\.\d+/); # Skip the header line
        my ($pos, $strand, $name) = (split /\t/, $i)[0,1,5];
        my ($start, $end) = (split /\.+/, $pos)[0,1];
        my $length = $end - $start + 1;
        $out{$name} = join "\t", ($chr, $length, $start, $end, $strand);
        $count ++;
    }
    close F2;
    pod2usage(-message => "input file not like a PTT file") if($count == 0);
    return %out;
}

# Read Sort
sub ReadSort{
    my $fh_sort = shift(@_);
    my %out = ();
    my $count = 0;
    open F, "<$fh_sort" or die "$!";
    while(my $i = <F>){
        chomp($i);
        my ($name, $chrom, $start, $end, $strand) = (split /\t/, $i)[0,1,3,4,5];
        next unless($start =~ /^\d+$/ || $end =~ /^\d+$/);
        my $length = $end - $start + 1;
        $out{$name} = join "\t", ($chrom, $length, $start, $end, $strand);
        $count ++;
    }
    pod2usage(-message => "Input file not like a Sort file") if($count == 0);
    return %out;
}

# Read FASTA
sub ReadFa{
    my $fh_fa = shift(@_);
    my $out = '';
    open F, "$fh_fa" or die "$!";
    my @lines = <F>;
    my $head = shift(@lines);
    my $id = (split /\s/,$head)[0];
    $id =~ s/\>//;
    my $seq = join '', @lines;
    $seq =~ s/\n//g;
    close F;
    return ($id,$seq);
}

# Txt2Seq
sub Txt2Seq{
    my ($seq, $start, $end, $strand) = @_;
    my $length = $end - $start + 1;
    $start -= 1; ## 1-st numbered as 0
    my $fa = substr($seq, $start, $length);
    $fa = &RevComp($fa) if($strand eq '-');

    my $fa_len = length($fa);
    my @seq_n = ();
    for(my $i = 0; $i<$fa_len; $i += 70){
        my $fa_n = substr($fa, $i, 70);
        push @seq_n, $fa_n;
    }
    my $fa_out = join "\n", @seq_n;
    return $fa_out;
}

# Reverse and complement seq
sub RevComp{
    my $in = shift(@_);
    my $out = reverse($in);
    $out =~ tr/ATCGatcg/TAGCtagc/;
    return $out;
}


=head1 NAME

c<sort2bed.pl> -- Transform the format between GFF/BED/PTT/Txt/FASTA... 

=head1 SYNOPSIS

perl sort2bed.pl -t sort2bed  -i input.bed -o out.bed

=head1 OPTIONS

=over 8

=item B<-t>, B<--trans>

Support the following transform:[gff|bed|sort|ptt]2fa, and 
[between 2 of these: gff,bed,sort,ptt]

=item B<-f>, B<--feature>

Select the feature that want to extract from GFF file: gene/exon/ncRNA/CDS/

=item B<-g>, B<--genome>

Choose the reference genome file

=item B<-o>, B<--output>

Write the result to the file. default:[inputname.]

=item B<-h>, B<--help>

Print this help

=item B<--man>

Print more detail help

=back

=head1 DESCRIPTION

This script was designed to perform format trnasformation between:
gff/bed/sort/ptt/txt/fasta

Explain these Format:
Find more details on UCSC website:  http://genome.ucsc.edu/FAQ/FAQformat.html

[GFF]    
    Have 9-required fields that must be tab-separated.

    1. seqname - The name of the seqeuence. Must be a chromosome or scaffold.       
    2. source - The program that generated this feature.
    3. feature - The name of this type of feature: CDS/start_codon/stop_codon/exon
    4. start - The starting position of the feature in sequence. The first base is numbered 1. *****
    5. end - The ending position of the feature.
    6. score - A score between 0 and 1000. use "." if there is no score.
    7. strand - "+" "-" or "." (don't know or don't care)
    8. frame - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. The value should be '.' if the feature is not a coding exon.
    9. group - All lines with the same group are linked together into a single item.

[BED]
    Have 3-required fileds and 9-additional optional fields

    1. chrom - The name of the chromosome or scaffold.
    2. chromStart - The starting position of the feature in the chromosome or scaffold. The first base is numbered 0.*****
    3. chromEnd - The ending position of the feature in the chromsome or scaffold.
    Optional fileds
    4. name - Define the name of the BED line.
    5. score - A score between 0 and 1000. Show the grey color.
    6. strand - Defines the strand. '+' or '-'.
    7. thickStart - The starting position at which the feature is drawn thickly (eg: start codon in gene). thickStart and thickEnd are usually set to chromStart when there is no thick part.
    8. thickEnd - eg: stop codon in gene.
    9. itemRgb - An RGB value of the form R,G,B (eg: 255,0,0). If the track line itemRgb is set ot "On".
    10. blockCount - The number of blocks (exons) in the BED line.
    11. blockSizes - A comma-separated list of the block sizes, correspond to blockCount.
    12. blockStarts - A comma-separated list of the block starts,

Example:
    browser position chr7:127471196-127495720
    browser hide all
    track name="ColorByStrandDemo" description="Color by strand demonstration"
    visibility=2 colorByStrand="255,0,0 0,0,255"
    chr7    127471196  127472363  Pos1  0  +
    chr7    127472363  127473530  Pos2  0  +
    chr7    127473530  127474697  Pos3  0  +
    chr7    127474697  127475864  Pos4  0  +
    chr7    127475864  127477031  Neg1  0  -
    chr7    127477031  127478198  Neg2  0  -
    chr7    127478198  127479365  Neg3  0  -
    chr7    127479365  127480532  Pos5  0  +
    chr7    127480532  127481699  Neg4  0  -

[PTT]
    The protein table file correspond with the CDS annotations from the GenBank file. Have 9 fileds.

    1. Location - The starting and ending position of the feature, separated by two "." (eg: 1..1024).
    2. Straind - "+" or "-".
    3. Length - The length of CDS.
    4. PID - The protein ID of the CDS.
    5. Gene - The gene name of the CDS, "-" for that have not a gene name.
    6. Synonym - The "locus_tag" from GFF file.
    7. Code - Not required.
    8. COG - The COG ID of the CDS.
    9. Product - The product of the CDS.

[Sorted]
    A tab-separated file with 6-required files and 6-additional option fields. This file store the information of ncRNAs.
    
    1. name - The name of the line.
    2. chrom - The name of the chromosome or scaffold.
    3. length - The length of the seq.
    4. start - The starting position of the seq. The first base is numbered 1. *****
    5. end - The ending position of the seq.
    6. strand - "+" or "-".
    7. pre_gene - The previous gene (left) of the seq.
    8. gap_1 - The distance between pre_gene and seq. (- to +)

[FASTA]
    fasta file.

=head1 AUTHOR

Wang Ming, wangmcas@gmail.com

=cut
