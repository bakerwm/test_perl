#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use File::Basename;
# Date: 19 Mar, 2013
# Update: 04 Jan, 2014
# Update: 14 Jan, 2014

sub help{
    print STDERR <<EOF
Usage: perl $0  input.gff  out.bed

Note:
    i       : The file in GFF3 format.
    o       : Output file in BED format. <input.bed>
EOF
}

die "Usage: perl $0 In.gff out.bed " unless(my $gff = shift);
my $BED = basename($gff); $BED =~ s/\.gff/.bed/;
my $outfile = "";
$outfile = $BED unless($outfile = shift);

open F,$gff or die "Cannot open file: $gff , $!\n";
my @gffs = <F>;
close F;

# Parsing the GFF file
my @Tabs = ();
foreach(@gffs){
    next unless(/\tgene\t/);
    my ($seqname, $start, $end, $strand) = (split /\t/)[0,3,4,6];
    my $len = $end - $start + 1;
    my ($GeneID, $GeneName, $Locus) = ();
    if(/Name=\w+/){ # NCBI gff version 1.20 & B42/BT annotation
        if(($GeneID, $Locus) = /GeneID\:(\d+).*locus_tag=(\w+)/){    # NCBI gff version 1.20
            ($GeneName) = /Name=(.*)\;Dbxref=/;
        }else{
            ($GeneName) = /Name=(\w+)\;/;
            $GeneID = $Locus = '-';  # B42/BT annotation
        }
    }else{
        ($Locus, $GeneID) = /locus_tag=(\w+)\;.*GeneID\:(\d+)/; # NCBI gff version 1.14
        if(($GeneName) = /\:(.*)\;locus_tag/){
        }else{
            $GeneName = '-';
        }
    }
   push @Tabs, (join"\t",($seqname, $start, $end, $GeneName, $len, $strand,$Locus));
}

open OUT,"> $outfile" or die "Cannot open file: $outfile, $!\n";
print OUT join"\n",(@Tabs),"\n";
close OUT;
