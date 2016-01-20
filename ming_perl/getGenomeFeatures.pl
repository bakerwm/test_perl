#!/usr/bin/perl -w
use strict;
use warnings;
use File::Path qw(make_path);
use File::Basename qw(basename dirname);
use File::Spec::Functions qw(catfile);
use Getopt::Std;

###########################################################################
# This script is designed to extract the following genomic features using #
# standard NCBI genome annotation files: GFF, PTT, RNT                    #
#                                                                         #
# The following features will report in BED format:                       #
# mRNA, asmRNA, tRNA, rRNA, IGR                                           #
#                                                                         #
# 2016-01-17                                                              #
###########################################################################

getGenomeFeatures();

sub getGenomeFeatures {
    my %opts = (s=>20);
    getopts("f:n:s:", \%opts);
    usage() if(@ARGV != 1);
    die("[-f] reference required\n") if( ! defined $opts{f} );
    my $igr_short = $opts{s};
    my $outdir = $ARGV[0];
    make_path($outdir) if( ! -d $outdir );
    ### parameters
    my $fna = $opts{f};
    my ($gff, $ptt, $rnt) = check_inputdir($fna);
    ### create reference list file
    my ($chr_name, $chr_length) = get_fnalist($fna);
    $chr_name = $opts{n} if(defined $opts{n});
    my $fna_list = catfile($outdir, "ref.list");
    open my $fh_list, "> $fna_list" or die "Cannot write to $fna_list, $!\n";
    print $fh_list $chr_name . "\t" . $chr_length . "\n";
    close $fh_list;
    ### get genomic features
    my $mRNA   = catfile($outdir, "mRNA.bed");
    my $asmRNA = catfile($outdir, "asmRNA.bed");
    my $tRNA   = catfile($outdir, "tRNA.bed");
    my $astRNA = catfile($outdir, "astRNA.bed");
    my $rRNA   = catfile($outdir, "rRNA.bed");
    my $asrRNA = catfile($outdir, "asrRNA.bed");
    my $igr    = catfile($outdir, "igr.bed");
    find_features($ptt, $rnt, $mRNA, $asmRNA, $tRNA, $astRNA, $rRNA, $asrRNA, $chr_name);
    ### igr
    extract_igr($mRNA, $tRNA, $rRNA, $igr, $outdir, $fna_list, $igr_short);
    print "Finish, get features!\n";
}

sub check_inputdir {
    my $fna = $_[0];
    my ($fna_name) = basename($fna) =~ /(.*)\.fna/;
    my $fna_dir    = dirname($fna);
    for my $feature ("fna", "gff", "ptt", "rnt") {
        my $check_file = catfile($fna_dir, $fna_name . "." . $feature);
        die("[$feature] file not found\n") if( ! -e $check_file );
    }
    my $gff = catfile($fna_dir, $fna_name . ".gff");
    my $ptt = catfile($fna_dir, $fna_name . ".ptt");
    my $rnt = catfile($fna_dir, $fna_name . ".rnt");
    return($gff, $ptt, $rnt);
}

sub get_fnalist {
    my $fna  = $_[0];
    my $chr_name;
    my $chr_length = 0;
    open my $fh_fna, "< $fna" or die "Cannot open fna file $fna, $!\n";
    while(<$fh_fna>) {
        chomp;
        if(/^\>/){
            $chr_name = (split /\s+/, $_)[0];
            if($chr_name =~ /\|/) {
                $chr_name = (split /\|/, $chr_name)[3];
            }
        }else {
            $chr_length += length($_);
        }
    }
    close $fh_fna;
    return($chr_name, $chr_length);
}

sub parse_ptt {
    my $ptt = $_[0];
    my %out = ();
    my $counter = 0;
    open my $fh_ptt, "< $ptt" or die "Cannot open $ptt, $!\n";
    while(<$fh_ptt>){
        chomp;
        next unless(/^\d+\.\.\d+/);
        my ($pos, $strand, $name) = (split /\t/)[0, 1, 5];
        my ($start, $end) = split /\.+/, $pos;
        my $length = $end - $start + 1;
        my $id = join(':', $name, $start, $end, $strand);
        $out{$id} = join("\t", $name, $length, $start, $end, $strand);
        $counter ++;
    }
    close $fh_ptt;
    die("No lines in ptt file\n") if($counter == 0);
    return(\%out);
}

sub parse_rnt {
    my $rnt      = $_[0];
    my %out = ();
    my $counter = 0;
    open my $fh_rnt, "< $rnt" or die "Cannot open $rnt, $!\n";
    while(<$fh_rnt>) {
        chomp;
        next unless(/^\d+\.\.\d+/);
        my ($pos, $strand, $name, $type) = (split /\t/)[0, 1, 5, 8];
        my ($start, $end) = split /\.+/, $pos;
        my $length = $end - $start + 1;
        my $id = $name . ":" . $type;
        $out{$id} = join("\t", $name, $length, $start, $end, $strand);
        $counter ++;
    }
    close $fh_rnt;
    die("No lines in rnt file\n") if($counter == 0);
    return(\%out);
}

sub find_features {
    my $ptt      = $_[0];
    my $rnt      = $_[1];
    my $mRNA     = $_[2];
    my $asmRNA   = $_[3];
    my $tRNA     = $_[4];
    my $astRNA   = $_[5];
    my $rRNA     = $_[6];
    my $asrRNA   = $_[7];
    my $chr_name = $_[8];
    ###
    my @mRNA_lists   = ();
    my @asmRNA_lists = ();
    my @tRNA_lists   = ();
    my @astRNA_lists = ();
    my @rRNA_lists   = ();
    my @asrRNA_lists = ();
    ### fetch mRNA and asmRNA
    my %for_ptt = %{ parse_ptt($ptt) };
    for my $p (keys %for_ptt) {
        my ($name, $length, $start, $end, $strand) = split /\t/, $for_ptt{$p};
        $start --;
        $end --;
        push @mRNA_lists, join("\t", $chr_name, $start, $end, $name, $length, $strand);
        my $as_name = 'AS_' . $name;
        my $as_strand = ($strand eq "+")?"-":"+";
        push @asmRNA_lists, join("\t", $chr_name, $start, $end, $as_name, $length, $as_strand);
    }
    ### fetch rRNA and tRNA, and asrRNA, astRNA
    my %for_rnt = %{ parse_rnt($rnt) };
    for my $r (keys %for_rnt) {
        my ($name, $length, $start, $end, $strand) = split /\t/, $for_rnt{$r};
        $start --;
        $end --;
        my $as_name = 'AS_' . $name;
        my $as_strand = ($strand eq "+")?"-":"+";
        if($r =~ /tRNA|Anticodon/) {
            push @tRNA_lists, join("\t", $chr_name, $start, $end, $name, $length, $strand);
            push @astRNA_lists, join("\t", $chr_name, $start, $end, $as_name, $length, $as_strand);
        }elsif($r =~ /ribosomal/) {
            push @rRNA_lists, join("\t", $chr_name, $start, $end, $name, $length, $strand);
            push @asrRNA_lists, join("\t", $chr_name, $start, $end, $as_name, $length, $as_strand);
        }else {
            next; ### !!!
        }
    }
    ### report
    open my $fh_mRNA,   "> $mRNA" or die "Cannot open $!\n";
    open my $fh_asmRNA, "> $asmRNA" or die "Cannot open $!\n";
    open my $fh_tRNA,   "> $tRNA" or die "Cannot open $!\n";
    open my $fh_astRNA, "> $astRNA" or die "Cannot open $!\n";
    open my $fh_rRNA,   "> $rRNA" or die "Cannot open $!\n";
    open my $fh_asrRNA, "> $asrRNA" or die "Cannot open $!\n";
    print $fh_mRNA   join("\n", @mRNA_lists) . "\n";
    print $fh_asmRNA join("\n", @asmRNA_lists) . "\n";
    print $fh_tRNA   join("\n", @tRNA_lists) . "\n";
    print $fh_astRNA join("\n", @astRNA_lists) . "\n";
    print $fh_rRNA   join("\n", @rRNA_lists) . "\n";
    print $fh_asrRNA join("\n", @asrRNA_lists) . "\n";
    close $fh_mRNA;
    close $fh_asmRNA;
    close $fh_tRNA;
    close $fh_astRNA;
    close $fh_rRNA;
    close $fh_asrRNA;
}

sub extract_igr {
    my $mRNA    = $_[0];
    my $tRNA    = $_[1];
    my $rRNA    = $_[2];
    my $igr     = $_[3];
    my $outdir  = $_[4];
    my $fa_list = $_[5];
    my $igr_short = $_[6];
    my $mtrRNA  = catfile($outdir, "mtrRNA.bed");
    my $igr_tmp = catfile($outdir, "igr_tmp.txt");
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
        next if($igr_length < $igr_short); ### skip the IGRs short than 20 bp
        my $id_p = sprintf"IGR%04d", 2 * $counter;
        my $id_n = sprintf"IGR%04d", 2 * $counter + 1;
        print $fh_igr join("\t", $chr, $start, $end, $id_p, $igr_length, "+") . "\n";
        print $fh_igr join("\t", $chr, $start, $end, $id_n, $igr_length, "-") . "\n";
        $counter ++;
    }
    close $fh_tmp;
    close $fh_igr;
    die("Find no igr fragements\n") if($counter == 0);
    unlink($igr_tmp) if( -e $igr_tmp );
    unlink($mtrRNA) if( -e $mtrRNA );
}

sub usage{
    die("
Usage: $0 [options] <outdir>

Options: -f  <STR>      : genome reference in fasta format, also require : 
                          [gff, ptt, rnt] files with the same name
         -n  <STR>      : the name of the genome
         -s  <INT>      : the smallest igr in length, default [20]

Note: 
    Only support 1-chromosome file; 
    Create IGR fragements based on mRNA/tRNA/rRNA

Example: 
    perl $0 -f NC_000962.fna output
\n");
}

###########################################################################
