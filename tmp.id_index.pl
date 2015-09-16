#!/usr/bin/env perl

###########################################################################
# Create id index for NCBI bactera genomes
#
# Input: GFF, PTT, FFN, FAA files
#
# Ouput: table contain id index for all genes
#
# Wang Ming wangmcas(AT)gmail.com
###########################################################################

use strict;
use warnings;
use File::Basename qw(basename dirname);
use File::Spec::Functions qw(catfile);
use File::Path qw(make_path);
use Getopt::Std;

use Data::Dumper;

id_index();
exit(1);

### subroutines
sub id_index {
    my %opts = ();
    getopts("o:", \%opts);
    usage() if(@ARGV != 1);
    my $in     = shift(@ARGV);
    my %ptt    = %{read_ptt($in)};
    die("[$in] not like a PTT file\n") if((keys %ptt) < 1 || ! $in =~ /\.ptt$/);
    my $strain   = `head -n1 $in`;
    $strain = (split /\-/, $strain)[0] if($strain =~ /\-/);
    ### parse other 3 files
    my ($g_name) = basename($in) =~ /(.*)\.ptt$/;
    my $outdir   = $g_name . "_index";
    $outdir      = $opts{o} if(defined $opts{o});
    make_path($outdir) if(! -d $outdir);
    my $GFF      = catfile(dirname($in), $g_name . ".gff");
    my $FFN      = catfile(dirname($in), $g_name . ".ffn");
    my $FAA      = catfile(dirname($in), $g_name . ".faa");
    die("[$in] file not exists\n") if(! -e $in);
    die("[*.gff, *.ffn, *.faa] missing some files\n") if(! -e $GFF || ! -e $FFN || ! -e $FAA);
    my %gff      = %{read_gff($GFF)};
    my %ffn      = %{read_ffn($FFN)};
    my %faa      = %{read_faa($FAA)};
    #
    my $gene_diff_pro = '';
    my $ffn_nothit    = '';
    my $faa_nothit    = '';
    my $report_gff  = 0;
    my $report_ptt  = 0;
    my $report_dna  = 0;
    my $report_pro  = 0;
    my $dna    = catfile($outdir, $g_name . "_dna.fa");
    my $pro    = catfile($outdir, $g_name . "_pro.fa");
    my $outidx = catfile($outdir, $g_name . "_index.txt");
    open my $fh_dna, "> $dna" or die "Cannot write to $dna, $!\n";
    open my $fh_pro, "> $pro" or die "Cannot write to $pro, $!\n";
    open my $fh_idx, "> $outidx" or die "Cannot write to $outidx, $!\n";
    for my $i (sort keys %{$ptt{index}}) {
        my ($id, $name, $length, $start, $end, $strand, $pid, $cog, $des) = split /\t/, $ptt{index}->{$i};
        ### search GFF
        my ($chr, $GeneID) = ('-', '-');
        if (exists $gff{index}->{$i}) {
            ($id, $chr, $length, $start, $end, $strand, $GeneID) = split /\t/, $gff{index}->{$i};
            delete $gff{index}->{$i};
            $report_gff ++;
        }elsif(exists $gff{locus}->{$id}) {
            ### output the records, gene and CDS have different coordinations
            my @tabs = split /\t/, $gff{locus}->{$id};
            $gene_diff_pro .= join("\t", @tabs[0..5], $id, $tabs[1],  $length, $start, $end, $strand, $des) . "\n";
            ($id, $chr, $length, $start, $end, $strand, $GeneID) = @tabs;
            $report_gff ++;
        }else {
            warn("[$id] not found in GFF file\n");
        }
        ### search FFN
        my $ffn_id = '-';
        my $ffn_fa = '';
        if(exists $ffn{$i}) {
            ($ffn_id, $ffn_fa) = split /\t/, $ffn{$i};
            print $fh_dna '>' . $id . "\n" . $ffn_fa . "\n";
            delete $ffn{$i};
            $report_dna ++;
        }else {
            $ffn_nothit .= join("\t", $id, $chr, $length, $start, $end, $strand, $des) . "\n";
#            die("[$id] not found in FFN file\n");
        }
        my $faa_id = '-';
        my $faa_fa = '';
        if(exists $faa{$pid}) {
            ($faa_id, $faa_fa) = split /\t/, $faa{$pid};
            print $fh_pro '>' . $id . "\n" . $faa_fa . "\n";
            delete $faa{$pid};
            $report_pro ++;
        }else {
             $faa_nothit .= join("\t", $id, $chr, $length, $start, $end, $strand, $des) . "\n";
#            die("[$id] not found in FAA file\n");
        }
        ### output : output the gene coordinations according to GFF file
        print $fh_idx join("\t", $id, $chr, $length, $start, $end, $strand, $name, $GeneID, 
                $pid, $cog, $ffn_id, $faa_id, $des) . "\n";
        $report_ptt ++;
    }
    close $fh_dna;
    close $fh_pro;
    close $fh_idx;
    ### validation
    
    ### report
    my $report_file = catfile($outdir, $g_name . '_report.txt');
    open my $fh_rpt, "> $report_file" or die "Cannot write to $report_file, $!\n";
    my $report_line = ">$g_name" . ': ' . $strain . "\n";
    $report_line   .= 'Find entries in PTT:' . "\t" . $report_ptt . "\n";
    $report_line   .= 'Find entries in FFN:' . "\t" . $report_dna . "\n";
    $report_line   .= 'Find entries in FAA:' . "\t" . $report_pro . "\n";
    $report_line   .= "\n" . 'The following gene(*.FFN) and CDS(*.FAA) are differ in coordinations:' . "\n";
    $report_line   .= join("\t", '#g.id', 'Chr', 'g.length', 'g.start', 'g.end', 'g.strand',
                           'p.id', 'Chr', 'p.length', 'p.start', 'p.end', 'p.strand', 'product') . "\n";
    $report_line   .= $gene_diff_pro . "\n";
    $report_line   .= "\n" . 'The following gene(s) are not found in FFN' . "\n";
    $report_line   .= $ffn_nothit . "\n";
    $report_line   .= "\n" . 'The following gene(s) are not found in FAA' . "\n";
    $report_line   .= $faa_nothit . "\n";
    print $report_line;
    print $fh_rpt $report_line;
    close $fh_rpt;
}

### read GFF
sub read_gff {
    my $in  = $_[0];
    my %gff = ();
    open my $fh_in, "< $in" or die "Cannot open $in, $!\n";
    while(<$fh_in>) {
        chomp;
        next if(/^\#|^\s*$/);
        next unless(/\tgene\t/);
        my $id       = '';
        my $GeneID   = '-';
        my @tabs     = split /\t/, $_;
        my $des      = $tabs[8];
        if(($id)     = $des =~ /locus_tag=(\w+\.*\w+)/){ ### WARN: multiple locus_tags 
        }elsif(($id) = $des =~ /Name=(\w+)/){
        }elsif(($id) = $des =~ /ID=gene\:(\w+)/){
        }else{
        }
        if(($GeneID) = $des =~ /GeneID\:(\d+)/){};
        # option
        my $chr      = $tabs[0];
        my $start    = $tabs[3];
        my $end      = $tabs[4];
        my $length   = $end - $start + 1;
        my $strand   = $tabs[6];
        my $index    = join("\:", $start, $end, $strand);
        my $out      = join("\t", $id, $chr, $length, $start, $end, $strand, $GeneID);
        $gff{index}->{$index} = $out;
        $gff{locus}->{$id}    = $out;
    }
    close $fh_in;
    return(\%gff);
}

### read PTT
sub read_ptt {
    my $in  = $_[0];
    my %ptt = ();
    open my $fh_in, "< $in" or die "Cannot open $in, $!\n";
    while(<$fh_in>) {
        chomp;
        next if(/^\#|^\s+$/);
        next unless(/^\d+\.\.\d+\t/);
        my @tabs = split /\t/, $_;
        my ($start, $end) = split /\.+/, $tabs[0];
        my $length = $end - $start + 1;
        my $strand = $tabs[1];
        my $pid    = $tabs[3];
        my $name   = $tabs[4];
        my $id     = $tabs[5];
        my $cog    = $tabs[7];
        my $des    = $tabs[8];
        my $index  = join("\:", $start, $end, $strand);
        $ptt{index}->{$index} = join("\t", $id, $name, $length, $start, $end, $strand, $pid, $cog, $des);
#        $ptt{pid}->{$pid} = join("\t", $id, $name, $length, $start, $end, $strand, $pid, $cog, $des);
    }
    close $fh_in;
    return(\%ptt);
}

### read FFN
sub read_ffn {
    my $in  = $_[0];
    my %ffn = ();
    open my $fh_in, "< $in" or die "Cannot open $in, $!\n";
    my @lines = <$fh_in>;
    close $fh_in;
    my @tabs  = split /\>/, join("", @lines);
    shift(@tabs);
    for my $i (@tabs) {
        my @aa = split /\n/, $i;
        my $header  = shift(@aa);
        my $fa      = join("\n", @aa);
        my $ffn_id  = (split /\s/, $header)[0];
        my ($a, $b) = $ffn_id =~ /\:(\w+)\-(\d+)/;
        my ($start, $end, $strand) = ($a, $b, '+');
        if($a =~ /^[a-z]/) {
            ($start) = $a =~ /(\d+)/;
            $end     = $b;
            $strand  = '-';
        }
        ($start, $end) = ($end, $start) if($start > $end);
        my $index      = join("\:", $start, $end, $strand);
        $ffn{$index}   = join("\t", $ffn_id, $fa);
    }
    return(\%ffn);
}

### read FAA
sub read_faa {
    my $in  = $_[0];
    my %faa = ();
    open my $fh_in, "< $in" or die "Cannot open $in, $!\n";
    my @lines = <$fh_in>;
    close $fh_in;
    my @tabs  = split /\>/, join("", @lines);
    shift(@tabs);
    for my $i (@tabs) {
        my @bb = split /\n/, $i;
        my $header = shift(@bb);
        my $fa     = join("\n", @bb);
        my $faa_id = (split /\s/, $header)[0];
        my ($pid, $np) = (split /\|/, $faa_id)[1,3];
        $faa{$pid}     = join("\t", $faa_id, $fa);
    }
    return(\%faa);
}

sub usage {
    die("
Usage: bac_index.pl [options] <in.ptt>

Options: -o <STR>       : save results to directory
         <in.ptt>       : contain GFF, FFN, FAA in the same directory

Output:
1. a script to retrieve id and sequence: input <GeneID, Locus, Name ...>, output <ID, sequence>

2. id_index.table
\n");
}

