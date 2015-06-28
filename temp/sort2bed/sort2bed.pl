#!/usr/bin/env perl
use strict;
use warnings;

##################################################
# This script is designed to transform the files
# between [GFF|BED|PTT|Sort|] format, and alos
# support read [blast tab] ouput (-m 8) and output
# Fasta file.
#
#
# 2015-05-01 wangmcas@gmail.com
##################################################

use File::Spec::Functions qw(catfile catdir);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path cwd);
use Getopt::Long;
use Pod::Usage;

my %opts = ();
my $type        = '';
my $input_file  = '';
my $output_file = '';
my $feature     = '';
my $genome      = '';
my $strain      = '';
my $help        = '';
my $man         = '';

GetOptions(
        't|type=s' => \$type,
        'i|input=s' => \$input_file,
        'o|output=s' => \$output_file,
        'f|feature=s' => \$feature,
        'g|genome=s' => \$genome,
        's|strain=s' => \$strain,
        'h|help' => \$help,
        'man' => \$man
        ) or pod2usage(2);

pod2usage(-verbose => 1) if($help);
pod2usage(-verbose => 2) if($man);
pod2usage(-message => "Need choose input file: [-i|--input]") if(! $input_file);

## Choose type of transformer
if(! $type){
    pod2usage(-message => "Need choose transform type: GFF|PTT|BED|Sort, eg: gff2bed");
}else{
    if($type =~ /^(gff|ptt|bed|sort|blast)2(gff|ptt|bed|sort|blast|fa)$/s){
        my ($t1, $t2) = $type =~ /^(\w+)2(\w+)$/;
        die "[-t $type] Choose another type\n" if ($t1 eq $t2);
    }else{
        pod2usage(-message => "Select [GFF,PTT,BED,Sort], as: gff2bed");
    }
}

## Check input file should be tab-spearated
&check_format($input_file);

##
my ($ta, $tb) = $type =~ /^(\w+)2(\w+)$/;

my %readin_data = ();
my %ourput_data = ();

# Original format
if($ta =~ /gff/s){
    pod2usage(-message => "Need input [-f|feature] ") if(! $feature);
    %readin_data = &read_gff($input_file, $feature);
}elsif($ta =~ /ptt/s){
    pod2usage(-message => "Need input [-g|genome] ") if(! $genome);
    %readin_data = &read_ptt($input_file, $genome);
}elsif($ta =~ /bed/s){
    %readin_data = &read_bed($input_file);
}elsif($ta =~ /sort/s){
    %readin_data = &read_sort($input_file);
}elsif($ta =~ /blast/s){
    %readin_data = &read_blast($input_file);
}else{
    # next;
}

$output_file = "out.txt" if(! $output_file);

open OUT, "> $output_file" or die "$!\n";
# Target format
if($tb =~ /gff/s){
    my @out = ();
    foreach my $n (sort keys %readin_data){
       my ($chr, $length, $start, $end, $strand) = split /\t/, $readin_data{$n};
       $feature = 'gene' if(! $feature);
       die "[$feature] not support: gene|exon|CDS|rRNA|tRNA|ncRNA\n" unless($feature =~ /gene|exon|CDS|rRNA|tRNA|ncRNA/);
       my $source = 'RefSeq';
       my $description = "ID=$n;gene=$n;Name=$n;gene_id=$n;locus_tag=$n";
       $chr = $strain if($strain);
       my $line = join "\t", ($chr, $source, $feature, $start, $end, '.', $strand, '.', $description);
       push @out, $line;
    }
    my $print_out = join "\n", @out;
    print OUT $print_out, "\n";
}elsif($tb =~ /ptt/s){
    my @out = ();
    foreach my $n (sort keys %readin_data){
        my ($chr, $length, $start, $end, $strand) = split /\t/, $readin_data{$n};
        $chr = $strain if($strain);
        my $colA = $start.'..'.$end;
        my $line = join "\t", ($colA, $strand, '-', '-', $n, $n, '-', '-', $chr);
        push @out, $line;
    }
    my $print_out = join "\n", @out;
    print OUT $print_out, "\n";
}elsif($tb =~ /bed/s){
    my @out = ();
    foreach my $n (sort keys %readin_data){
        my ($chr, $length, $start, $end, $strand, $tail) = split /\t/, $readin_data{$n},6;
#        $tail = '-' unless($tail);
## 1-st base numbered 0
        $start -= 1;
        $end -= 1;
        $chr = $strain if($strain);
        my $line = join "\t", ($chr, $start, $end, $n, $length, $strand, $tail);
        push @out, $line;
    }
    my $print_out = join "\n", @out;
    print OUT $print_out, "\n";
}elsif($tb =~ /sort/s){
    my @out = ();
    foreach my $n (sort keys %readin_data){
        my ($chr, $length, $start, $end, $strand, $tail) = split /\t/, $readin_data{$n},6;
#        $tail = '-' unless($tail);
        $chr = $strain if($strain);
        push @out, join "\t", ($n, $chr, $length, $start, $end, $strand, $tail);
    }
    my $print_out = join "\n", @out;
    print OUT $print_out, "\n";
}elsif($tb =~ /fa/s){
    pod2usage(-message => "Need reference file [-g|genome]") if(! $genome);
    my ($id, $ref) = &read_fasta($genome);
    my @out = ();
    foreach my $n (sort keys %readin_data){
        my ($chr, $length, $start, $end, $strand) = split /\t/, $readin_data{$n};
        my $seq = &txt_to_seq($ref, $start, $end, $strand);
        push @out, (">$n", $seq);
    }
    my $print_out = join "\n", @out;
    print OUT $print_out, "\n";
}else{
}

### subroutines ###
## Check tab separated file [gff/ptt/BED/sort]
sub check_format{
    my $f = shift(@_);
    open my $f_chk, "< $f" or die "$!\n";
    while (my $i = <$f_chk>){
        next if($i =~ /^\#/);
        if($i =~ /\t/){ # at least 6 columns
            die "Input is not tab-separated: $!\n" if((split /\t/,$i) < 6);
        }
    }
    close $f_chk;
}

## Readin GFF
sub read_gff {
    my ($fh_gff, $feature) = @_;
    my %tabs = ();
    my $count = 0;
    open my $gff, "<$fh_gff" or die "$!\n";
    while(my $g = <$gff>){
        next if($g =~ /^\s*$/); # skip blank lines
        next if($g =~ /^\#/);
        chomp ($g);
        next unless ($g =~ /\t$feature\t/);
        my ($chrom, $ref, $fe, $start, $end, $t1, $strand, $t2, $des, $tail) = split /\t/, $g,10;
        my $length = $end - $start + 1;
        my $name = '-';
        $tail = '' unless($tail);
        # feature = gene
        if($feature eq 'gene') {
            if(($name) = $des =~ /locus_tag=(\w+)/){
            }elsif(($name) = $des =~ /Name=(\w+)/){
            }elsif(($name) = $des =~ /ID=gene\:(\w+)/){
            }elsif(($name) = $des =~ /GeneID\:(\d+)/){
            }else{
            }
        # feature = CDS, exon, rRNA, tRNA, ncRNA
        }elsif($feature =~ /CDS|rRNA|tRNA|ncRNA|exon/){
            if(($name) = $des =~ /gene=(\w+)/) {
            }elsif(($name) = $des =~ /ID=\w+\:(\w+)/ ){
            }elsif(($name) = $des =~ /GeneID\:(\d+)/){
            }elsif(($name) = $des =~ /Name=(\w+\-*\d)\;/){
            }else{
            }
        # feature not found
        }else{
            pod2usage(-message => "[$feature] is not supported.");
        }
        $tabs{$name} = join "\t", ($chrom, $length, $start, $end, $strand, $tail);
        $count ++;
    }
    close $gff;
    pod2usage(-message => "[$feature] not found in GFF file") if($count == 0);
    return %tabs;
}

# Read BED
sub read_bed{
    my $fh_bed = shift(@_);
    my %out = ();
    my $count = 0;
    open F, "<$fh_bed" or die "$!\n";
    while(my $i = <F>){
        next if ($i =~ /^\s*$/); # skip blank lines
        next if($i =~ /^\#/);
        chomp($i);
        next if($i =~ /^\s*$/); # Skip blank lines
        my ($chrom, $start, $end, $name, $cont, $strand, $tail) = split /\t/, $i, 7;
        $tail = ' ' unless($tail);
        $start += 1; ## BED first base = 0;
        $end += 1; ## BED first base = 0;
        my $length = $end - $start + 1;
        $out{$name} = join "\t", ($chrom, $length, $start, $end, $strand, $tail);
        $count ++;
    }
    close F;
    pod2usage(-message => "[$input_file] not a BED file") if($count == 0);
    return %out;
}

# Read PTT
sub read_ptt{
    my ($fh_ptt, $g) = @_;
    my %out = ();
    my $count = 0;
    my $chr = '';
    open F, $g or die "Cannot open $g, $!\n";
    my $line = <F>;
    $chr = (split /\s/, $line)[0]; $chr =~ s/\>//;
#    ($chr) = (split /\|/, $line)[3];
    close F;
    open F2, "<$fh_ptt" or die "$!\n";
    while(my $i = <F2>){
        next if ($i =~ /^\s*$/); # skip blank lines
        next if($i =~ /^\#/);
        chomp($i);
        next unless($i =~ /^\d+\.\.\d+/); # Skip the header line
        my ($pos, $strand, $name) = (split /\t/, $i)[0,1,5];
        my ($start, $end) = (split /\.+/, $pos)[0,1];
        my $length = $end - $start + 1;
        $out{$name} = join "\t", ($chr, $length, $start, $end, $strand, '');
        $count ++;
    }
    close F2;
    pod2usage(-message => "[$input_file] not a PTT file") if($count == 0);
    return %out;
}

# Read Sort
sub read_sort{
    my $fh_sort = shift(@_);
    my %out = ();
    my $count = 0;
    open F, "<$fh_sort" or die "$!\n";
    while(my $i = <F>){
        next if ($i =~ /^\s*$/); # skip blank lines
        next if($i =~ /^\#/);
        chomp($i);
        $i =~ s/\r//; # chomp the Windown newline
        my ($name, $chrom, $cont, $start, $end, $strand, $tail) = split(/\s+/, $i, 7);
        $tail = ' ' unless($tail);
#        my ($name, $chrom, $start, $end, $strand) = (split /\t/, $i)[0,1,3,4,5];
        next unless($start =~ /^\d+$/ || $end =~ /^\d+$/);
        my $length = $end - $start + 1;
        $out{$name} = join("\t", $chrom, $length, $start, $end, $strand, $tail);
        $count ++;
    }
    pod2usage(-message => "[$input_file] not a Sort file") if($count == 0);
    return %out;
}

## Readin blast
sub read_blast{
    my $fh_blast = shift(@_);
    my %out = ();
    my $count = 0;
    open F, "< $fh_blast" or die "$!\n";
    while(my $i = <F>){
         next if ($i =~ /^\s*$/); # skip blank lines
         next if($i =~ /^\#/);
         chomp($i);
         my ($q_name, $chrom, $s_start ,$s_end) = (split /\t/, $i)[0,1,8,9];
         my $strand = '+';
         ($s_start, $s_end, $strand) = ($s_end, $s_start, '-') if ($s_start > $s_end);
         my $length = $s_end - $s_start + 1;
         $out{$q_name} = join "\t", ($chrom, $length, $s_start, $s_end, $strand, '');
         $count ++;
    }
    pod2usage(-message => "[$input_file] not a blast -m 8 output") if($count == 0);
    return %out;
}

# Read FASTA
sub read_fasta{
    my $fh_fa = shift(@_);
    my $out = '';
    open F, "$fh_fa" or die "$!\n";
    my @lines = <F>;
    my $head = shift(@lines);
    my $id = (split /\s/,$head)[0];
    $id =~ s/\>//;
    my $seq = join '', @lines; $seq =~ s/\n//g;
    close F;
    return ($id,$seq);
}

# txt_to_seq
sub txt_to_seq{
    my ($seq, $start, $end, $strand) = @_;
    my $length = $end - $start + 1;
    $start -= 1; ## 1-st numbered as 0
    my $fa = substr($seq, $start, $length);
    $fa = &rev_comp($fa) if($strand eq '-');
    return $fa;
}

# Reverse and complement seq
sub rev_comp{
    my $in = shift(@_);
    my $out = reverse($in);
    $out =~ tr/ATCGatcg/TAGCtagc/;
    return $out;
}

=head1 NAME

c<sort2bed.pl> -- Transform the format between GFF/BED/PTT/Sort/Blast/FASTA...

=head1 SYNOPSIS

perl sort2bed.pl -t sort2bed  -i input.bed -o out.bed

=head1 OPTIONS

=over 8

=item B<-t>, B<--trans>

Support the following transform:[gff|bed|sort|ptt|blast]2fa, and
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
    A tab-separated file with 6-required fields and 6-additional option fields. This file store the information of ncRNAs.

    1. name - The name of the line.
    2. chrom - The name of the chromosome or scaffold.
    3. length - The length of the seq.
    4. start - The starting position of the seq. The first base is numbered 1. *****
    5. end - The ending position of the seq.
    6. strand - "+" or "-".
    7. pre_gene - The previous gene (left) of the seq.
    8. gap_1 - The distance between pre_gene and seq. (- to +)
    9. nex_gene - The next gene (right) of the seq.
    10. gap_2 - The distance between seq and nex_gene.
    11. direction - the direction of the three seq. eg: /+/-/+/
    12. description - the description of the seq. AS = antisense, IGR = intergenic, IM = within mRNA, PM = part of a mRNA.

[Blast]
    A tab-separated file with at least 12-required fields. The -m 8 output of blastn.

    1. q_name   - The query name.
    2. s_name   - The subject name, also the chromosome name.
    3. Identity - The identity of the hit.
    4. Length   - The lenght of the hit.
    5. mismatch - The number of mismatches in the alignment.
    6. gap      - The number of gaps in the alignment.
    7. q_start  - The start of the query.
    8. q_end    - The end of the query.
    9. s_start  - The start of the subject.
    10. s_end   - The end of the subject.
    11. E-value - The e-value of the hit.
    12. score   - The bit score of the hit.

[FASTA]
    fasta file.

=head1 AUTHOR

Wang Ming, wangmcas@gmail.com

=cut


################
# Change log
# 2015-6-19: delete windows "return line", s/\r//; in read_sort
