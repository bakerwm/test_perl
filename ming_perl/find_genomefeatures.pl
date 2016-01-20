#!/usr/bin/perl -w
use strict;
use warnings;
use File::Path qw(make_path);
use File::Spec::Functions qw(catfile);
use Getopt::Std;

###########################################################################
# This script is designed to extract the following genomic features using #
# standard NCBI genome annotation files: GFF, PTT, RNT                    #
#                                                                         #
# The following features will report in BED format:                       #
# mRNA, asmRNA, tRNA, rRNA, IGR                                           #
#                                                                         #
###########################################################################

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
my $mtrRNA;
my $gene;
my $chr_name;

parse_para();
find_features();

sub parse_para{
    my %opts = ();
    getopts("a:n:", \%opts);
    usage() if(@ARGV != 1);
    $out_path = $ARGV[0];
    make_path($out_path) if( ! -d $out_path);
    die("[-a] reference file required\n") if( ! defined $opts{a} );
    my $fa_name = $fa = $opts{a};
    $fa_name =~ s/\.f[n]*a//;
    $gff = $fa_name . '.gff';
    $ptt = $fa_name . '.ptt';
    $rnt = $fa_name . '.rnt';
    foreach my $f ("gff", "ptt", "rnt"){
        my $ff = $fa_name . "." . $f;
        die("[$f] file not found") if( ! -e $ff);
    }
    $mRNA = catfile($out_path, "mRNA.bed");
    $as   = catfile($out_path, "asmRNA.bed");
    $igr  = catfile($out_path, "IGR.bed");
    $tRNA = catfile($out_path, "tRNA.bed");
    $rRNA = catfile($out_path, "rRNA.bed");
    $mtrRNA = catfile($out_path, "mtrRNA.bed");
    $gene = catfile($out_path, "gene.bed");
    ### parse reference file
    $fa_list = catfile($out_path, "ref.txt");
    my $chr_length = 0;
    ($chr_name, $chr_length) = read_fa();
    $chr_name = $opts{n} if defined $opts{n};
    open my $fh_list, "> $fa_list" or die "Cannot write to fa_list file, $fa_list, $!\n";
    print $fh_list $chr_name . "\t". $chr_length . "\n";
    close $fh_list;
}

sub find_features{
    open my $fh_mRNA, "> $mRNA" or die "$!";
    open my $fh_as,   "> $as"   or die "$!";
    open my $fh_igr,  "> $igr"  or die "$!";
    open my $fh_tRNA, "> $tRNA" or die "$!";
    open my $fh_rRNA, "> $rRNA" or die "$!";
    # mRNA & AS
    my %p = read_ptt();
    foreach my $i (keys %p){
        my ($chr, $length, $start, $end, $strand, $name) = split /\t/,$p{$i};
        $start--; # bed file 0-left most
        $end--;
        print $fh_mRNA join("\t", ($chr_name, $start, $end, $name, $length, $strand)), "\n";
        my $as_str = ($strand eq '+')?'-':'+';
        print $fh_as join("\t", ($chr_name, $start, $end, "AS_".$name, $length, $as_str)), "\n";
    }
    close $fh_mRNA;
    close $fh_as;
    # rRNA & tRNA
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
    # IGR
    extract_igr();
}

sub extract_igr {
    my $igr_tmp = catfile($out_path, "igr_tmp.txt");
    system("cat $mRNA $tRNA $rRNA | sort -k1,1 -k2,2n > $mtrRNA");
    system("bedtools complement -i $mtrRNA -g $fa_list > $igr_tmp");
    open my $fh_tmp, "< $igr_tmp" or die "Cannot open $igr_tmp, $!\n";
    open my $fh_igr, "> $igr" or die "Cannot open $igr, $!\n";
    my $counter = 0;
    while(<$fh_tmp>) {
        chomp;
        my ($chr, $start, $end) = split /\t/, $_;
        $start ++;
        $end --;
        my $igr_length = $end - $start + 1;
        next if($igr_length < 20); ### skip the IGRs short than 20 bp
        my $id_p = sprintf"IGR%04d", 2 * $counter;
        my $id_n = sprintf"IGR%04d", 2 * $counter + 1;
        print $fh_igr join("\t", $chr, $start, $end, $id_p, $igr_length, "+") . "\n";
        print $fh_igr join("\t", $chr, $start, $end, $id_n, $igr_length, "-") . "\n";
        $counter ++;
    }
    close $fh_tmp;
    close $fh_igr;
    unlink($igr_tmp) if( -e $igr_tmp );
    unlink($mtrRNA) if( -e $mtrRNA );
}

sub read_fa{
    my $chr_name;
    my $chr_length;
    my $header;
    open my $fh_fa, "< $fa" or die "$!";
    while(<$fh_fa>){
        chomp;
        if(/^\>/){
            $_ =~ s/^\>//;
            $header = $_;
        }else{
            $chr_length += length($_); 
        }
    }
    close $fh_fa;
    $chr_name = (split /\s+/, $header)[0];
    if($header =~ /\|/) {
        $chr_name = (split /\|/, $header)[3];
    }
    return ($chr_name, $chr_length);
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
        my $id = join("\:", ($chr_name, $start, $end, $strand));
        $out{$id} = join "\t", ($chr_name, $length, $start, $end, $strand, $name);
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
        my $line = join("\t", ($chr_name, $start, $end, $name, $length, $strand));
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

