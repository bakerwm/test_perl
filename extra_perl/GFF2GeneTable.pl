#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use File::Basename;
# Date: 19 Mar, 2013
# Upudate: Line 49, $GeneName for Rvnt
#
#
sub help{
    print STDERR <<EOF
Usage: perl $0  -i  input.gff  -o  out.table

Note:
    i       : The file in GFF3 format.
    o       : Output file with 7 lines in each row.
              <GeneID> <Seqname> <Start> <End> <Strand> <GeneName> <Locus>
              default: input.table
EOF
}

my %options = ();
getopt("i:o:",\%options);
if(!defined $options{i} ){
    &help();
    exit(1);
}

my $filename = basename($options{i});
my $outname  = $filename;
$outname =~ s/\.gff/\.table/; # outfile name default is input.table

$options{o} = (defined $options{o})?$options{o}:"$outname";
my $Table = '';
open F,$options{i} or die "Cannot open file: $options{i} , $!\n";
while(<F>){
    chomp;
    next unless(/\tgene\t/);
    my @ps = split(/\t/, $_);
    my $seqname = $ps[0];
    my ($start, $end, $strand) = ($ps[3], $ps[4], $ps[6]);
    my ($GeneName, $GeneID, $Locus) = qw(- - -);
    ($GeneName) = ($ps[8]=~/locus_tag=(\w+)\;/);

    if($ps[8] =~ /Name.*gbkey.*locus_tag/){
## !gff-spec-version 1.20
        $_ =~ /Name=(.*)\;Dbxref=GeneID\:(\d+)\;.*locus_tag=(\w+)/;
        ($GeneName, $GeneID, $Locus) = ($1, $2, $3);
        ($GeneName) = ($GeneName =~ /note=(\w+\-\w+)/)?$1:$GeneName;  # update: 2014-11-04
#        if($ps[8]=~/Rvnt/){
#            $ps[8] =~ /note=(\w+\-\w+)/i;
#            $GeneName = $1;
#        }

    }else{
## !gff-spec-version 1.14
        $ps[8] =~ /locus_tag=(\w+)\;.*GeneID:(\d+)/;
        ($Locus, $GeneID) = ($1, $2);
        $ps[8] =~ /ID=\w+.\d\:(\w+)\;locus_tag/;
        $GeneName = ($1)?$1:'-';
        if($ps[8]=~/Rvnt/){
            $ps[8] =~ /note=(\w+\-\w+)/;
            $GeneName = $1;
        }
    }
#    my @gids = split(/\;/,$GeneID) 
    if($GeneName =~ /\;/){
        my @t = split(/\;/,$GeneName);
        $GeneName = $t[0];
    }

    my $out = join"\t",($GeneID, $seqname, $start, $end, $strand, $GeneName, $Locus);
    $Table .= $out."\n";
}
close F;

my @head = ('GeneID', 'seqname', 'start', 'end', 'strand', 'GeneName', 'Locus');
my $headOut = join"\t",(@head);
$Table = $headOut."\n".$Table;

open OUT,"> $options{o}" or die "Cannot open file: $options{o}, $!\n";
print OUT $Table;
