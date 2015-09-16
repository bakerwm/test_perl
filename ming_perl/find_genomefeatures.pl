#!/usr/bin/perl -w
use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec::Functions qw(catfile);
use Getopt::Std;

my $fa;
my $fa_list;
my $gff;
my $ptt;
my $rnt;
my $out_path;
my $mRNA;
my $as;
my $igr;
my $rRNA;
my $tRNA;
my $gene;
my $chrom;

parse_para();
find_features();
#exit(1);

sub parse_para{
    my %opts = ();
    getopts("a:n:", \%opts);
#    my $usage = "Usage: $0 [-a] <ref.fa> [-n] <ref.name> out.dir\n";
#    die("$usage") if(@ARGV != 1);
    usage() if(@ARGV != 1);
    die("[-a] need ref fa\n") unless defined $opts{a};
    my $fa_name = $fa = $opts{a};
    $fa_name =~ s/\.f[n]*a//;
    $gff = $fa_name . '.gff';
    $ptt = $fa_name . '.ptt';
    $rnt = $fa_name . '.rnt';
    foreach my $f ($fa, $gff, $ptt, $rnt){
        die("[$f] file not found") unless -e $f;
    }
    $out_path = shift(@ARGV);
    make_path("$out_path") unless -d $out_path;    
    $mRNA = catfile($out_path, "mRNA.bed");
    $as   = catfile($out_path, "AS.bed");
    $igr  = catfile($out_path, "IGR.bed");
    $tRNA = catfile($out_path, "tRNA.bed");
    $rRNA = catfile($out_path, "rRNA.bed");
    $gene = catfile($out_path, "gene.bed");
    $fa_list = catfile($out_path, "ref.txt");
    
    my ($header, $length) = read_fa();
    $chrom = $header;
    $chrom = $opts{n} if defined $opts{n};
    open my $fh_list, "> $fa_list" or die "$!";
    print $fh_list $chrom . "\t". $length . "\n";
    close $fh_list;
}

sub find_features{
    open my $fh_mRNA, "> $mRNA" or die "$!";
    open my $fh_as, "> $as" or die "$!";
    open my $fh_igr, "> $igr" or die "$!";
    open my $fh_tRNA, "> $tRNA" or die "$!";
    open my $fh_rRNA, "> $rRNA" or die "$!";
    open my $fh_gene, "> $gene" or die "$!";
    # mRNA & AS
    my %p = read_ptt();
    foreach my $i (keys %p){
        my ($chr, $length, $start, $end, $strand, $name) = split /\t/,$p{$i};
        $start--; # bed file 0-left most
        $end--;
        print $fh_mRNA join("\t", ($chrom, $start, $end, $name, $length, $strand)), "\n";
        my $as_str = ($strand eq '+')?'-':'+';
        print $fh_as join("\t", ($chrom, $start, $end, "AS_".$name, $length, $as_str)), "\n";
    }
    close $fh_mRNA;
    close $fh_as;
    # IGR
    my %g = read_gff();
    foreach my $n (keys %g){
        my ($chr, $length, $start, $end, $strand, $name) = split /\t/, $g{$n};
        $start --;
        $end --;
        print $fh_gene join("\t", ($chrom, $start, $end, $name, $length, $strand)), "\n";
    }
    close $fh_gene;
    system "bedtools complement -i $gene -g $fa_list > tmp.bed";
    # filter IGR file, and trim 1-base at both ends
    system "cat tmp.bed | awk \'{\$2++; \$3--; if((\$3-\$2)>20) print \$1\"\\t\"\$2\"\\t\"\$3 }\' > tmp_trim.bed  ";
    open my $fh_tmp, "< tmp_trim.bed" or die "$!";
    my $counter = 1;
    while(<$fh_tmp>){
        chomp;
        my ($chr, $start, $end) = split /\t/;
        my $length = $end - $start;
        my $id_p = sprintf"IGR%04d", $counter;
        $counter ++;
        my $id_n = sprintf"IGR%04d", $counter;
        $counter ++;
        print $fh_igr join("\t", ($chrom, $start, $end, $id_p, $length, '+')), "\n";
        print $fh_igr join("\t", ($chrom, $start, $end, $id_n, $length, '-')), "\n";
    }
    close $fh_tmp;
    close $fh_igr;
    unlink("tmp.bed", "tmp_trim.bed", $gene);
    # rRNA & tRNa
    my %t = read_rnt();
    foreach my $m (keys %t){
        if($m =~ /tRNA|Anticodon/){
            print $fh_tRNA $t{$m}, "\n";
        }elsif($m =~ /ribosomal\sRNA/){
            print $fh_rRNA $t{$m}, "\n";
        }else{
        }
    }
    close $fh_tRNA;
    close $fh_rRNA;
}

sub read_fa{
    my $length;
    my $header;
    open my $fh_fa, "< $fa" or die "$!";
    while(<$fh_fa>){
        chomp;
        if(/^\>/){
            $_ =~ s/^\>//;
            $header = (split /\s/, $_)[0];
        }else{
            $length += length($_); 
        }
    }
    close $fh_fa;
    return ($header, $length);
#    $chrom = $opts{n} if defined $otps{n};
#    open my $fh_list, "> $fa_list" or die "$!";
#    print $fh_list join("\t", ($chrom, $length)), "\n";
#    close $fh_list;
}

sub read_gff{
    my %tabs = ();
    my %ncRNA = ();
    my $count = 0;
    open my $fh_gff, "< $gff" or die "$!";
    while(my $g = <$fh_gff>){
        chomp ($g);
        next if($g =~ /^\#/);
        next unless($g =~ /\t(gene|ncRNA)/);
        my ($chr, $start, $end, $strand, $des) = (split /\t/, $g)[0,3,4,6,8];
        my $length = $end - $start + 1;
        my $name = '-';
        if(($name) = $des =~ /locus_tag=(\w+)/){
        }elsif(($name) = $des =~ /Name=(\w+)/){
        }elsif(($name) = $des =~ /ID=gene\:(\w+)/){
        }elsif(($name) = $des =~ /GeneID\:(\d+)/){
        }else{
        }
        my $id = join("\:", ($chr,$start, $end, $strand));
        my $line = join "\t", ($chr, $length, $start, $end, $strand, $name);
        if($g =~ /\tgene\t/){
            $tabs{$id} = $line;
        }elsif($g =~ /\tncRNA/){
            $ncRNA{$id} = $line;
        }
        $count ++;
    }
    close $fh_gff;
    die("not found gene in GFF file") if($count == 0);
    foreach my $i (keys %ncRNA){
        delete $tabs{$i};
    }
    return %tabs;
}

sub read_ptt{
    my %out;
    my $count = 0;
    open my $fh_ptt, "< $ptt" or die "$!";
    while(my $p = <$fh_ptt>){
        chomp ($p);
        next unless($p =~ /^\d+\.\.\d+/); # Skip the header line
        my ($pos, $strand, $name) = (split /\t/, $p)[0,1,5];
        my ($start, $end) = (split /\.+/, $pos)[0,1];
        my $length = $end - $start + 1;
        my $id = join("\:", ($chrom, $start, $end, $strand));
        $out{$id} = join "\t", ($chrom, $length, $start, $end, $strand, $name);
        $count ++;
    }
    close $fh_ptt;
    die("found no lines in PTT file") if($count == 0);
    return %out;
}

sub read_rnt{
    my %out;
    my $count = 0;
    open my $fh_rnt, "< $rnt" or die "$!";
    while(my $p = <$fh_rnt>){
        chomp ($p);
        next unless($p =~ /^\d+\.\.\d+/); # Skip the header line
        my ($pos, $strand, $name, $type) = (split /\t/, $p)[0,1,5,8];
        my ($start, $end) = (split /\.+/, $pos)[0,1];
        my $length = $end - $start + 1;
        $start --;
        $end --;
        my $line = join("\t", ($chrom, $start, $end, $name, $length, $strand));
        my $id = $name . ':' . $type;
        $out{$id} = $line;
        $count ++;
    }
    close $fh_rnt;
    die("found no lines in PTT file") if($count == 0);
    return %out;
}

sub usage {
    die("
Usage: find_genomefeatures.pl [options] <outdir>

Options: -a : <STR>     : fasta file of NCBI bacteria genome
                          need GFF, PTT, RNT files in the same dir
         -n : <STR>     : name of Chr for ouput
Note: 
1. only output IGRs > 20 nt

Example: 
find_genomefeatures.pl -a ~/NCBI/NC_000962.fna -n NC_000962.3 out
\n");
}

