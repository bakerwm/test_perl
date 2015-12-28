#!/usr/bin/env perl

############################################################
# Add annotations from multiple files to a text.
#
# Wang Ming wangmcas(AT)gmail.com
# 2015-07-01
##############################################################

use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;

addNotes();
exit (1);

###
sub addNotes {
    my %opts = (n => 1);
    getopts("i:n:l:o:", \%opts);
    usage() if(@ARGV == 0);
    die("[-i $opts{i}] input file not exists\n") if(! -e $opts{i});
    my %info = %{ parseNotes(@ARGV) };
    my $col_n = $opts{n} - 1; # perl 0-lest most index
    open my $fh_in, "< $opts{i}" or die "Cannot open $opts{i}, $!\n";
    my $fh_out = *STDOUT;
    if(defined $opts{o}) {
        open $fh_out, "> $opts{o}" or die "Cannot open $opts{o}, $!\n";
    }
    my $total  = 0;
    my $update = 0;
    while(<$fh_in>) {
        chomp;
        my @tabs = split(/\s+/, $_);
        if(! $opts{n} =~ /^\d+$/ || $opts{n} > @tabs) {
            die("[col:$opts{n}] out range of [$opts{i}]\n");
        }
        my ($tail, $dd) = combineNotes($tabs[$col_n], \%info);
        print $fh_out join("\t", @tabs, $tail) . "\n";
        if($dd) {
            $update ++;
        }
        $total ++;
    }
    if(defined $opts{l}) {
        open my $fh_l, "> $opts{l} " or die "Cannot write to $opts{l}, $!\n";
        print $fh_l "Update: $update of $total\n";
        close $fh_l
    }
    close $fh_in;
    close $fh_out;
}

sub combineNotes {
    my $id = $_[0];
    my %hf = %{$_[1]};
    my @tags = ();
    my $count = 0;
    for my $i (sort keys %hf) {
        my $flag = '-';
        if(exists $hf{$i}->{$id}) {
            $flag = join("\,", @{$hf{$i}->{$id}});
            $count ++
        }
        push @tags, $flag;
    }
    my $n = join("\t", @tags);
    return ($n, $count);
}

sub parseNotes {
    my @lists = @_; # @ARGV
    my %info  = ();
    my $count = 1;  # file order
    my $file  = '';
    my $id    = 1;
    my $note  = 2;
    for my $m (@lists) {
        my ($file, $a, $b1, $b2) = split /\t/, check_fmt($m);
        # perl index (0-start)
        $b1 --;
        $b2 --;
        open my $fh_in, "< $file" or die "Cannot open $file, $!\n";
        while(<$fh_in>) {
            chomp;
            my @tabs = split(/\s+/, $_);
            my $col_id = $a - 1;
            if($a <= @tabs && $b1 <= @tabs && $b2 <= @tabs) {
                push @{$info{$count}->{$tabs[$col_id]}}, join("\t", @tabs[$b1..$b2]);
            }else {
                die("[$a,$b1,$b2] the column numbers are out of range\n");
            }
        }
        close $fh_in;
        # make sure: duplicates not allowed when add "multiple columns" to file
        if(($b2 - $b1) > 0) {
            for my $n (keys %{$info{$count}}) {
                if( @{$info{$count}->{$n}} > 1 ) {
                    die("[$n, $b1\-$b2] multiple hits found, not allowed when multiple columns selected to add\n");
                }
            }
        }
        $count ++;
    }
    return(\%info);
}

# format: file,m,1-3
sub check_fmt {
    my $in = $_[0];
    my $flag = 0;
    if( (split/\,/, $in) == 3 ) {
        my($file, $a, $b) = split /\,/, $in;
        die("[$file] file not exists\n") if(! -e $file);
        my $b1 = my $b2 = $b;
        ($b1, $b2) = split /\-/, $b if($b =~ /\d+\-\d+/);
        if($a=~/^\d+$/ && $b1=~/^\d+$/ && $b2=~/^\d+$/ && $b1 <= $b2){
           $flag = join("\t", $file, $a, $b1, $b2);
        }
    }
    if( ! $flag ) {
        die("[$in] is expected in : <file>,<N>,<N1>-<N2>\n");
    }
    return $flag;
}

sub usage {
    die(qq/
Usage: addNotes.pl [options] file1,1 file2,2 file3,2 ...

Add Notes from other files according to the id filed.

Options: -i     : The input txt file [id in column-1]
         -n     : The column contain the id field. [1]
         -o     : (optional) if specified, redirect output to file
                  [STDOUT]

Format: <FILE>,<N1>,<N2>

FILE: the input file
N1   : the column of id in FILE
N2   : columns of info to add, <INT> or <INT-INT>, to specify the note info in FILE

Example:
1. add info in col-7 of fileA to f.txt
perl addNotes.pl -i f.txt -o f.info fileA,1,7

2. add info in col-7 to col-10 of fileA to f.txt
perl addNotes.pl -i f.txt -o f.info fileA,1,7-10
\n/);
}

__END__

### Change log

2015-07-14
  1. add: "one query have multi hits in info"
  2. add option: "-l" redirect log to file
