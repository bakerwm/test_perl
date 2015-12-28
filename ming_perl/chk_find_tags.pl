#!/usr/bin/env perl

############################################################################
# This script is designed to pick out the regions with continous expression 
# that covered with reads more than a cutoff.
# 1. edges > cutoff
# 2. length > 20
# 3. edge cov no les than 1/2 of mean_cov #
# 4. edge no less that 1/10 of peak cov #
#
# Author: Wang Ming (wangmcas@gmail.com)
# April 21, 2015
############################################################################

use warnings;
use strict;
use Getopt::Std;
use Data::Dumper;

sub usage {
    die("Usage: find_tags.pl [options] <cov.file>

Options: -c STR         : set a cutoff to search the tags, default [100]
         -s STR         : the strand of iput file, (+/-), default [+]
\n");
}

my %opts = (c => 100, s => '+');
getopts("c:s:", \%opts);
usage() if(@ARGV != 1);

my $strand = $opts{"s"};
my $cutoff = $opts{"c"};
my $cutoff_length = 20;
my $flag = '_p' if ($strand eq '+');
$flag = '_n' if ($strand eq '-');

my %cov = ();
my %cord = ();
my $genome_length = 1;
my $chr_name = '';

my $infile = $ARGV[0];
open my $fh, "<$infile" or die "$!";
while(<$fh>){
    chomp;
    my ($chr, $pos, $count) = split /\t/;
    $chr_name = $chr;
    $cord{$pos} = $count;
    $cov{$count}->{$pos} = 1 if ($count >= $cutoff);
    $genome_length = ($genome_length < $pos)?$pos:$genome_length;
}
close $fh;

my %hited = ();
my $counter = 1;
foreach my $m (sort {$b<=>$a} keys %cov) {
    foreach my $p (sort {$a<=>$b} keys %{$cov{$m}} ) {
        next if (exists $hited{$p});
        my ($left, $left_cov, $right, $right_cov) =  &search_ends($p);
        # check errors
        if($left eq '' || $right eq '') {
        }
        # skip empty lines
        next unless(defined $left_cov && defined $right_cov);
        my $id_er = sprintf "%04d", $counter;
        my $seq_length = $right - $left + 1;
        my $seq_line = join "\t", ("Seed$id_er$flag", $chr_name, $seq_length, $left,
                                   $right, $strand, $left_cov, $m, $right_cov);
        # skip short seq
        next if($seq_length < $cutoff_length);
        # mean cov to edges
        my $mean_cov = &mean_cov($left, $right);
        next if($mean_cov < 2 * $cutoff);
        # edge_cov > 100 * flank_cov
        print $seq_line, "\t", $mean_cov, "\n";
        $counter ++;
    }
}

### Subroutines ###
sub mean_cov {
    my ($left, $right) = @_;
    my $length = $right - $left + 1;
    my $sum_cov = 0;
    for(my $i = $left; $i <= $right; $i++) {
        next unless(exists $cord{$i});
        $sum_cov += $cord{$i};
    }
    my $mean_cov = sprintf "%.4f", $sum_cov/$length;
    return $mean_cov;
}

sub search_ends {
    my $pos  = $_[0];
    my $left = $pos;
    my $left_cov = $cord{$pos};
    my $right = $pos;
    my $right_cov = $cord{$pos};
    my @leftinfos = ($left, $left_cov);
    my @rightinfos = ($right, $right_cov);
    # search left
    for(my $i = $pos; $i >= 1; $i--) {
        my $i_cov = $cord{$i};
        my $i_cov_next = $cord{$i + 1};
        my $i_cov_pre = $cord{$i - 1};
        $i_cov_next = 1 if($i_cov_next == 0);
        $i_cov_pre = 1 if($i_cov_pre == 0);
        @leftinfos = ($i, $i_cov);
        if(exists $hited{$i - 1} || $i == 1){
            last;
        }else{
            if($i_cov >= $cutoff){
                if(($i_cov/$i_cov_pre) > 100){ # folder max/edge
                    last;
                }else{
                    $hited{$i} = 1;
                }
            }else{
                @leftinfos = (($i + 1), $i_cov_next);
                last;
            }
        }
    }
    # search right
    for(my $n = $pos; $n < $genome_length; $n++) {
        my $n_cov =$cord{$n};
        my $n_cov_next = $cord{$n + 1};
        my $n_cov_pre = $cord{$n - 1};
        $n_cov_next = 1 if($n_cov_next == 0);
        $n_cov_pre = 1 if($n_cov_pre == 0);
        @rightinfos = ($n, $n_cov);
        if(exists $hited{$n + 1} || $n == $genome_length) {
            last; 
        }else{
            if($n_cov >= $cutoff) {
                if(($n_cov/$n_cov_next) > 100) { # folder max/edge
                    last;
                }else{
                    $hited{$n} = 1;
                }
            }else{
                @rightinfos = (($n - 1), $n_cov_pre);
                last;
            }
        }
    }
    return (@leftinfos, @rightinfos);
}

#### STOP
__END__

=head1 NAME

c<get_sRNA_cutoff.pl> - Find sRNA transcrpts according to the coverage file.

=head1 SYNOPSIS

perl get_sRNA_cutoff.pl -i Rv.cov.n -c 100 -s + -o Rv.temp

=head1 OPTIONS

=over 8

=item B<-i> str, B<--input>

The tab-separated file contain the reads coverage depth at each base of the genome.
<Strain/Chr> <Position> <Count>

=item B<-c> N, B<-cutoff>

The base covered with at least cut-off reads will be checked for sRNA transcripts.
recommend [100]

=item B<-s> str, B<--strand>

The strand of the input file +/-. default [+]

=item B<-o> str, b<--output>

The file contain the sRNA transcripts. default [out.txt]

=item B<-h>, B<--help>

Print this help.

=head1 AUTHOR

Wang Ming (wangmcas@gmail.com)

