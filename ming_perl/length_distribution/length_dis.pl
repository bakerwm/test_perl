#!/usr/bn/perl

##--------------------------------
## length_dis.pl - statistic the length distribution of input file
##
## Ming Wang wangmcas@gmail.com, 2016-06-01
## 
## This program is distributed in the hope that it will be usefull, but 
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABLITY 
## or FITNESS. 
## 
## FOR A PARTICULAR PURPOSE. See the GNU General Public License version 3 
## for more details.
##--------------------------------

use strict;
use warnings;
use File::Basename qw(basename);
use Cwd;
use Getopt::Std;

use Data::Dumper;

# Usage
sub usage {
    print "
length_dis.pl - statistic the length distribution of input file

Usage: perl $0 [options]

  required:
    -i      input file, support: fastq, fasta, BED, tab-separated file ... 
    -t      type of input file: [fasta, fastq, bed, tab, ...]
    -o      output txt file
  
  optional:
    -v
\n";
    exit 0;
}

# option
my %opts = ();
getopts("i:o:t:c:v:h", \%opts);

usage if(defined $opts{h});
die "[-i] need to specify input file (-h for more details)\n" if(! defined $opts{i});
die "[-t] need to specify type of input file (-h for more details)\n" if(! defined $opts{t});
die "[-o] need to spedify the output file\n" if(! defined $opts{o});
my $infile  = $opts{i};
my $outfile = $opts{o}; #basename($infile) . ".LenDis.txt";
my $col     = $opts{c}; # for TAB input, choose the length where
my $count_in_id = (defined $opts{c})?$opts{c}:0; # for BED input, count in id?

# main script
if ($opts{t} =~ /^fastq|fq$/i) {
    awk_fq($infile, $outfile);
}elsif ($opts{t} =~ /^fasta|fa$/i) {
    awk_fa($infile, $outfile);
}elsif ($opts{t} =~ /^bed$/i) {
    parse_bed($infile, $outfile, $count_in_id); 
}elsif ($opts{t} =~ /^tab$/i) {
    die "[-c] need to be specified, if using TAB input\n" if(! defined $opts{c});
    die "[-c] have to be Integer\n" if(! $opts{c} =~ /^\d+$/);
    parse_tab($infile, $outfile, $col);
}else{
    die "[$opts{t}] unknown type of input, see -h for more help\n";
}

# subroutines
sub awk_fq {
    my $infile  = $_[0];
    my $outfile = $_[1];
    # FastQ
    `awk 'NR%4 == 2 {lengths[length(\$0)]++} END {for (l in lengths) {print l, lengths[l]}}' $infile > $outfile `
}

sub awk_fa {
    my $infile  = $_[0];
    my $outfile = $_[1];
    # FastQ
    `awk 'NR%2 == 0 {lengths[length(\$0)]++} END {for (l in lengths) {print l, lengths[l]}}' $infile > $outfile `
}

sub parse_bed {
    my $infile  = $_[0];
    my $outfile = $_[1];
    my $count_in_id = $_[2];
    #
    open my $f_in, "< $infile" or die "Cannot open $infile, $!\n";
    open my $f_out, "> $outfile" or die "Cannot open $outfile, $!\n";
    my %ha = ();
    while(<$f_in>) {
        chomp;
        my @tabs = split ("\t", $_);
        next if(@tabs < 6); ## at least 6 columns
        next if(! $tabs[1] =~ /^\d+$/); ## col-2, start site
        next if(! $tabs[2] =~ /^\d+$/); ## col-3, end site
        my $length = $tabs[2] - $tabs[1]; ## length
        my $count = 1;
        if($count_in_id) { ## read count within id
            my $count = (split "\,", $tabs[3])[-1];
            $count =~ s/\s+//;
            $count = ($count =~ /^\d+$/)?$count:1;
        }
        $ha{$tabs[5]}->{$length} += $count;
    }
    #
    for my $strand (keys %ha) {
        for my $length (sort {$a<=>$b} keys %{$ha{$strand}}) {
            print $f_out join ("\t", $length, $ha{$strand}->{$length}, $strand) . "\n";
        }
    }
    close $f_in;
    close $f_out;
}

sub parse_tab {
    my $infile   = $_[0];
    my $outfile  = $_[1];
    my $col      = $_[2]; # which col contain, seq or length
    #
    $col --; ## 0-start in perl
    open my $f_in, "< $infile" or die "Cannot open $infile, $!\n";
    open my $f_out, "> $outfile" or die "Cannot opne $outfile, $!\n";
    my %ha = ();
    while(<$f_in>) {
        chomp;
        my @tabs = split ("\t", $_);
        if(@tabs < $col) {
            print STDERR "Warning: line-$$ contains fields less than -c\n";
        }
        ##
        my $length = 0;
        my $hit = $tabs[$col];
        if ($hit =~ /^(A|T|C|G)*$/) { # sequence
            $length = length($hit);
        }elsif ($hit =~ /^\d+$/) {
            $length = $hit;
        }else {
            print STDERR "Warning: line-$$ contains neither sequence nor Integer in column-[$col]\n";
        }
        $ha{$length} ++;
    }
    #
    for my $len (sort {$a<=>$b} keys %ha) {

        next if($len == 0);
        print $f_out join("\t", $len, $ha{$len}) . "\n";
    }
    close $f_in;
    close $f_out;
}

# EOF #
