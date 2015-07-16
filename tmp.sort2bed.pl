#!/usr/bin/env perl

##################################################
# This script is designed to transform the files
# between [GFF|BED|PTT|Sort|] format, and alos
# support read [blast tab] ouput (-m 8) and output
# Fasta file.
#
# 2015-05-01 wangmcas@gmail.com
# update: 2015-07-16
##################################################

use strict;
use warnings;
use Getopt::Long;

sort2bed();
exit(1);

###
sub sort2bed {
    my %opts = ();
    my $from;
    my $to;
    my $outfile;
    my $feature_type = 'gene';
    my $ref;
    my $chr    = '';
    my $extra  = 1;
    my $help;
    my $man;
    GetOptions(
            'a|fmtin=s'  => \$from,
            'b|fmtout=s' => \$to,
            'o|output=s' => \$outfile,
            't|type=s'   => \$feature_type,
            'f|ref=s'    => \$ref,
            'n|chr=s'    => \$chr,
            'x|extra=i'  => \$extra,
            'h|help'     => \$help,
            'man'        => \$man,
            ) or usage_smp();
    usage_smp()  if($help);
    usage_full() if($man);
    usage_smp()  if(@ARGV == 0 && -t STDIN);
    die("[-a, -b] both are required: GFF,BED,PTT,Sort,Fasta\n") if(! defined $from || ! defined $to);
    check_type($from, $to);
    ###
    my $fh_out = *STDOUT;
    if(defined $opts{o}) {
        open $fh_out, "> $opts{o}" or die "Cannot write to $opts{o}, $!\n";
    }
    $chr = (defined $chr)?$chr:0;
    while(<>) {
        chomp;
        next if(/^\#|^\s*$/);
        next if((split /\t/) < 6);
        my $fmt_in = guess_fmt($_);
        if($fmt_in ne $from) {
            die("[-a $from] format error, is it <$fmt_in>?\n$_\n");
        }
        my $out = fmt_convert($_, $from, $to, $feature_type, $ref, $chr, $extra);
        next if(! $out);
        print $fh_out $out . "\n";
    }
}
 
sub check_type {
    my $in  = $_[0];
    my $out = $_[1];
    if( ! $in =~ /gff|ptt|bed|sort|blast/ ) { # && ! $out =~ /gff|ptt|bed|sort|blast|fa/ ) {
        die("[-a $in] unknown format\n");
    }elsif( ! $out =~ /gff|ptt|bed|sort|blast|fa/ ) {
        die("[-b $out] unknown format\n");
    }
}

### guess the format of input
sub guess_fmt {
# Sort, BED PTT, GFF, Blast    
    my $in  = $_[0];
    my $fmt = '';
    my @tabs = split /\t/, $in;
    if(@tabs >= 6) {
        if($in =~ /\d+\t\d+\t[+-]/) {
            $fmt = 'sort';
        }elsif($in =~ /\d+\t\d+\t.*\t\d+\t[+-]/) {
            $fmt = 'bed';
        }elsif($in =~ /^\d+\.\.\d+\t/) {
            $fmt = 'ptt';        
        }elsif($in =~ /\d+\t\d+\t\.\t[+-]\t\.\t/) {
            $fmt = 'gff';
        }elsif($in =~ /\d+\.(\d+\t){7}/) {
            $fmt = 'blast8';
        }else {
            $fmt = 'na';
        }
    }else {
        $fmt = 'na';
    }
    return $fmt;
}

#############################
# Convert format
#############################
sub fmt_convert {
    my $line = $_[0];
    my $from = $_[1];
    my $to   = $_[2];
    my $type = $_[3];
    my $ref  = $_[4];
    my $name = $_[5];
    my $k    = $_[6]; # control tail
    my $in;
    my $out;
    ###
    if($from =~ /^gff$/i) {
        $in = read_gff($line, $type, $name);
    }elsif($from =~ /^ptt$/i) {
        die("[-f] need input the reference fasta file\n") if(! defined $ref);
        $in = read_ptt($line, $name, $ref);
    }elsif($from =~ /^bed$/i) {
        $in = read_bed($line, $name);
    }elsif($from =~ /^sort$/) {
        $in = read_sort($line, $name);
    }elsif($from =~ /^blast$/) {
        $in = read_blast8($line, $name);
    }else {
        die("[-a $from] unknown format\n");
    }
###
    if($to =~ /^gff$/i) {
        $out = write_gff($in, $type, $k);
    }elsif($to =~ /^bed$/) {
        $out = write_bed($in, $k);
    }elsif($to =~ /^sort$/) {
        $out = write_sort($in, $k);
    }elsif($to =~ /^ptt$/) {
        $out = write_ptt($in, $k);
    }elsif($to =~ /^fa$/) {
        die("[-f] Need input reference file, fasta\n") if(! defined $ref);
        $out = write_fa($in, $ref, $k);
    }else {
        die("[-b $to] unknown format\n");
    }
    return $out; 
}

### read GFF
sub read_gff {
    my ($line, $type, $name) = @_;
    my $out  = '';
    my $tail = '';
    my $gene_id = '';
    my @tabs = split /\t/, $line, 10;
    $tail    = $tabs[9] if(defined $tabs[9]);
    my $des  = $tabs[8];
    if($type eq $tabs[2]) {
        if($type eq 'gene') {
            if(($gene_id) = $des =~ /locus_tag=(\w+)/){
            }elsif(($gene_id) = $des =~ /Name=(\w+)/){
            }elsif(($gene_id) = $des =~ /ID=gene\:(\w+)/){
            }elsif(($gene_id) = $des =~ /GeneID\:(\d+)/){
            }else{
            }
        # feature = CDS, exon, rRNA, tRNA, ncRNA
        }elsif($type =~ /CDS|rRNA|tRNA|ncRNA|exon/){
            if(($gene_id) = $des =~ /gene=(\w+)/) {
            }elsif(($gene_id) = $des =~ /ID=\w+\:(\w+)/ ){
            }elsif(($gene_id) = $des =~ /GeneID\:(\d+)/){
            }elsif(($gene_id) = $des =~ /Name=(\w+\-*\d)\;/){
            }else{
            }
        # feature not found
        }else{
            die("[-t $type] unknown type\n");
        }
        # option
        my $chr = ($name)?$name:$tabs[0];
        my $start  = $tabs[3];
        my $end    = $tabs[4];
        my $length = $end - $start + 1;
        my $strand = $tabs[6];
        $out = join("\t", $gene_id, $chr, $length, $start, $end, $strand, $tail);
    }
    return $out;
}

### read BED
sub read_bed {
    my $line    = $_[0];
    my $name    = $_[1]; # chr name
    my $out     = '';
    my $tail    = '';
    my @tabs    = split /\t/, $line, 7;
    $tail       = $tabs[6] if(defined $tabs[6]);
    my $start   = $tabs[1];
    my $end     = $tabs[2];
    my $length  = $end - $start + 1;
    my $gene_id = $tabs[3];
    my $strand  = $tabs[5];
    my $chr     = ($name)?$name:$tabs[0];
    # bed is 0-left most index
    $start ++;
    $end ++;
    $out = join("\t", $gene_id, $chr, $length, $start, $end, $strand, $tail);
    return $out;
}

### read PTT
sub read_ptt {
    my $line = $_[0];
    my $name = $_[1]; # chr name
    my $ref  = $_[2];
    my $ref_id = parse_refid($ref);
    my $out  = '';
    my $tail = '';
    my @tabs = split /\t/, $line;
    if($line =~ /^\d+\.\.\d+/) {
        my $tail    = $tabs[4] . "\t" . $tabs[8];
        my ($start, $end) = split /\.+/, $tabs[0];
        my $length  = $end - $start + 1;
        my $strand  = $tabs[1];
        my $gene_id = $tabs[5];
        my $chr     = ($name)?$name:$ref_id;
        $out = join("\t", $gene_id, $chr, $length, $start, $end, $strand, $tail);
    }
    return $out;
}

sub parse_refid {
    my $in = $_[0];
    my $header = qx(grep '>' $in);
    my $id;
    if ($header =~ /\s/) {
        my $h1 = (split /\s/, $header)[0];
        if($h1 =~ /\|/) {
            $id = (split /\|/, $h1)[3];
        }
    }
    return $id;
}

# read sort
sub read_sort {
    my $line = $_[0];
    my $name = $_[1];
    my $out  = '';
    my $tail = '';
    my @tabs = split /\t/, $line, 7;
    $tail = $tabs[6] if(defined $tabs[6]);
    my $gene_id = $tabs[0];
    my $chr    = ($name)?$name:$tabs[1];
    my $start  = $tabs[3];
    my $end    = $tabs[4];
    my $strand = $tabs[5];
    my $length = $end - $start + 1;
    $out = join("\t", $gene_id, $chr, $length, $start, $end, $strand, $tail);
    return $out;
}

### read blast8 
sub read_blast8 {
    my $line = $_[0];
    my $name = $_[1];
    my $out  = '';
    my $tail = '';
    my @tabs = split /\t/, $line;
    my ($q_name, $subject, $s_start, $s_end) = (split /\t/, $line)[0,1,8,9];
    my $strand;
    ($s_start, $s_end, $strand) = ($s_start < $s_end)?($s_start, $s_end, '+'):($s_end, $s_start , '-');
    my $length = $s_end - $s_start + 1;
    my $chr = ($name)?$name:$subject;
    $out = join("\t", $q_name, $chr, $length, $s_start, $s_end, $strand, $tail);
    return $out;
}

### read fasta file, for 1-chr strain
sub read_fasta {
    my $in  = $_[0];
    my $out = '';
    open my $fh_in, "< $in" or die "Cannot open $in, $!\n";
    my @lines = <$fh_in>;
    my $id = (split /\s/, shift(@lines))[0];
    $id =~ s/\>//;
    my $seq = join'', @lines;
    $seq =~ s/\n//g;
    close $fh_in;
    return $seq;
}

sub txt_to_seq {
    my $in     = $_[0];
    my $start  = $_[1];
    my $end    = $_[2];
    my $strand = $_[3];
    $start --; # perl index 0-left most
    my $length = $end - $start + 1;
    my $fa = substr($in, $start, $length);
    if($strand eq '-') {
        $fa =~ tr/ATCGagcg/TAGCtagc/;
        $fa = reverse($fa);
    }
    return $fa;
}

sub write_gff { # not include tail in GFf output
    my $in = $_[0];
    my $type = $_[1];
    my $k  = $_[2]; # control output tail
    my $source = 'RefSeq';
    my $out = '';
    if($in) {
        my ($id, $chr, $length, $start, $end, $strand, $tail) = split /\t/, $in,7;
        my $description = "ID=$id;Name=$id;gene_id=$id;locus_tag=$id";    
        $out  = join("\t", $chr, $source, $type, $start, $end, '.', $strand, '.', $description);
    }
    return $out;
}

sub write_bed {
    my $in = $_[0];
    my $k  = $_[1];
    my $out = '';
    if($in) {
        my ($id, $chr, $length, $start, $end, $strand, $tail) = split /\t/, $in,7;
        $start --;
        $end --;
        $out = join("\t", $chr, $start, $end, $id, $length, $strand);
        $out .= "\t" . $tail if($k); # 
    }
    return $out;
}

sub write_ptt {
    my $in = $_[0];
    my $k  = $_[1];
    my $out = '';
    if($in) {
        my ($id, $chr, $length, $start, $end, $strand, $tail) = split /\t/, $in,7;
        my $p  = $start . '..' . $end;
        $out = join("\t", $p, $strand, '-', '-', $id, $id, '-', '-', $chr);
        $out .= "\t" . $tail if($k);
    }
    return $out;
}

sub write_sort {
    my $in = $_[0];
    my $k  = $_[1];
    my $out = '';
    if($in) {
        my ($id, $chr, $length, $start, $end, $strand, $tail) = split /\t/, $in,7;
        $out = join("\t", $id, $chr, $length, $start, $end, $strand);
        $out .= "\t" . $tail if($k);
    }
    return $out;
}

sub write_fa {
    my $in  = $_[0];
    my $ref = $_[1];
    my $k   = $_[2];
    my $out = '';
    if($in) {
        my ($id, $chr, $length, $start, $end, $strand, $tail) = split /\t/, $in,7;
        my $ref_fa = read_fasta($ref);
        my $fa = txt_to_seq($ref_fa, $start ,$end, $strand);
        $out = '>' . $id . "\n" . $fa;
    }
    return $out;
}

sub usage_smp {
    die("
Usage: tab_convert.pl [options] <input|STDIN>

Converter the format of input between (GFF, BED, PTT, Sort, Blast, Fasta);

Opeiont: -o --output    <STR>       : if specified, redirect results to file. (optional)
         -a --fmtin     <STR>       : the format of input file: GFF, BED, PTT, Sort, Blast, [sort]
         -b --fmtout    <STR>       : the format of output results, GFF, BED, PTT, Sort, Fasta, [bed]
         -t --type      <STR>       : the type of feature in GFF file. [gene]
         -f --ref       <STR>       : the reference in fasta format.
         -n --chr       <STR>       : the name of chr.
         -x --extra     <INT>       : whether report the tail in result. 0=no, 1=yes [1];
         -h --help                  : show this help
         -man                       : show more details

<input|STDIN> : have to be Tab-separated.

Example:
tab_convert.pl -a gff -b sort -t gene NC_123456.gff > strain.txt
tab_convert.pl -a ptt -b sort -f ref.fa NC_123456.ptt > strain.txt
tab_convert.pl -n NC_123456 -x 0 in.txt > out.bed # exclude the tail in output
\n");
}

sub usage_full {
    die("
Usage: tab_convert.pl [options] <input|STDIN>

Converter the format of input between (GFF, BED, PTT, Sort, Blast, Fasta);

Opeiont: -o --output    <STR>       : if specified, redirect results to file. (optional)
         -a --fmtin     <STR>       : the format of input file: GFF, BED, PTT, Sort, Blast, [sort]
         -b --fmtout    <STR>       : the format of output results, GFF, BED, PTT, Sort, Fasta, [bed]
         -t --type      <STR>       : the type of feature in GFF file. [gene]
         -f --ref       <STR>       : the reference in fasta format.
         -n --chr       <STR>       : the name of chr.
         -x --extra     <INT>       : whether report the tail in result. 0=no, 1=yes [1];
         -h --help                  : show this help
         -man                       : show more details

<input|STDIN> : have to be Tab-separated.

Example:
tab_convert.pl -a gff -b sort -t gene NC_123456.gff > strain.txt
tab_convert.pl -a ptt -b sort -f ref.fa NC_123456.ptt > strain.txt
tab_convert.pl -n NC_123456 -x 0 in.txt > out.bed # exclude the tail in output

Description:
The program was designed to convert the file between the following formats:
GFF, BED, PTT, Sort, Fasta, Blast

Find more details about file format on UCSC website: http://genome.ucsc.edu/FAQ/FAQformat.html

[GFF]
    Have 9-required fields that must be tab-separated.

    1. seqname - The name of the seqeuence. Must be a chromosome or scaffold.
    2. source  - The program that generated this feature.
    3. feature - The name of this type of feature: CDS/start_codon/stop_codon/exon
    4. start   - The starting position of the feature in sequence. The first base is numbered 1. *****
    5. end     - The ending position of the feature.
    6. score   - A score between 0 and 1000. use "." if there is no score.
    7. strand  - "+" "-" or "." (don't know or don't care)
    8. frame   - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. The value should be '.' if the feature is not a coding exon.
    9. group   - All lines with the same group are linked together into a single item.

[BED]
    Have 3-required fileds and 9-additional optional fields

    1. chrom        - The name of the chromosome or scaffold.
    2. chromStart   - The starting position of the feature in the chromosome or scaffold. The first base is numbered 0.*****
    3. chromEnd     - The ending position of the feature in the chromsome or scaffold.
    Optional fileds
    4. name         - Define the name of the BED line.
    5. score        - A score between 0 and 1000. Show the grey color.
    6. strand       - Defines the strand. '+' or '-'.
    7. thickStart   - The starting position at which the feature is drawn thickly (eg: start codon in gene). thickStart and thickEnd are usually set to chromStart when there is no thick part.
    8. thickEnd     - eg: stop codon in gene.
    9. itemRgb      - An RGB value of the form R,G,B (eg: 255,0,0). If the track line itemRgb is set ot \"On\".
    10. blockCount  - The number of blocks (exons) in the BED line.
    11. blockSizes  - A comma-separated list of the block sizes, correspond to blockCount.
    12. blockStarts - A comma-separated list of the block starts,

Example:
    browser position chr7:127471196-127495720
    browser hide all
    track name=\"ColorByStrandDemo\" description=\"Color by strand demonstration\"
    visibility=2 colorByStrand=\"255,0,0 0,0,255\"
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

    1. Location - The starting and ending position of the feature, separated by two \".\" (eg: 1..1024).
    2. Straind  - \"+\" or \"-\".
    3. Length   - The length of CDS.
    4. PID      - The protein ID of the CDS.
    5. Gene     - The gene name of the CDS, \"-\" for that have not a gene name.
    6. Synonym  - The \"locus_tag\" from GFF file.
    7. Code     - Not required.
    8. COG      - The COG ID of the CDS.
    9. Product  - The product of the CDS.

### this is a customised format for small RNA analysis
[Sorted] 
    A tab-separated file with 6-required fields and 6-additional option fields. This file store the information of ncRNAs.

    1. name         - The name of the line.
    2. chrom        - The name of the chromosome or scaffold.
    3. length       - The length of the seq.
    4. start        - The starting position of the seq. The first base is numbered 1. *****
    5. end          - The ending position of the seq.
    6. strand       - \"+\" or \"-\".
    7. pre_gene     - The previous gene (left) of the seq.
    8. gap_1        - The distance between pre_gene and seq. (- to +)
    9. nex_gene     - The next gene (right) of the seq.
    10. gap_2       - The distance between seq and nex_gene.
    11. direction   - the direction of the three seq. eg: /+/-/+/
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

Wang Ming, wangmcas(AT)gmail.com
2015-07-16
\n");
}

__END__

Change log
2010-4-15
  v0.1 Start this program

2012-4-6
  v0.1 Add more features, support GFF, BED, Sort formats

...

2014-11-5
  v0.5 support GFF, BED, Sort, fasta, blast formats

2015-6-19
  v1.0 delete windows "return line", s/\r//; in read_sort

2015-07-16
  v1.1 rewrite the program, support STDIN and STDOUT
