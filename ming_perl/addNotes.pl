#!/usr/bin/env perl

############################################################
# Add annotations from multiple files to a text.
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
    getopts("i:n:o:", \%opts);
    usage() if(@ARGV == 0);
    die("[-i $opts{i}] input file not exists\n") if(! -e $opts{i});
    die("[-o] Need input file of output\n") if(! defined $opts{o});

    my $i = parseNotes(@ARGV);
    my %info = %{$i};
    my $col_n = $opts{n} - 1; # perl 0-lest most index
    open my $fh_in, "< $opts{i}" or die "Cannot open $opts{i}, $!\n";
    open my $fh_out, "> $opts{o}" or die "Cannot open $opts{o}, $!\n";
    while(<$fh_in>) {
        chomp;
        my @tabs = split(/\s+/, $_);
        my $tail = combineNotes($tabs[$col_n], \%info);
        print $fh_out join("\t", @tabs, $tail) . "\n";
    }
    close $fh_in;
    close $fh_out;
}

sub combineNotes {
    my ($id, $h) = (@_);
    my %hf = %{$h};
    my @tags = ();
    for my $i (sort keys %hf) {
        my $flag = '-';
        if(exists $hf{$i}->{$id}) {
            $flag = $hf{$i}->{$id};
        }
        push @tags, $flag;
    }
    return join("\t", @tags);
}

sub parseNotes {
    my @lists = @_;
    my %info  = ();
    my $count = 1;
    my $file  = '';
    my $id    = 1;
    my $note  = 2;
    for my $m (@lists) {
        die("[$m] format not match: [<file>,<num>,<num>]\n") if(! check_fmt($m) );
        my ($file, $a, $b) = check_fmt($m);
#print STDERR join("\n", $m, $file, $a, $b) . "\n";
        my ($b1, $b2) = $b =~ /(\d+)\:(\d+)/;
        # perl index (0-start)
        $b1 --;
        $b2 --;
        open my $fh_in, "< $file" or die "Cannot open $file, $!\n";
        while(<$fh_in>) {
            chomp;
            my @tabs = split(/\s+/, $_);
            my $col_id = $a - 1;
            $info{$count}->{$tabs[$col_id]} = join("\t", @tabs[$b1..$b2]);
        }
        close $fh_in;
        $count ++;
    }
    return(\%info);
}

# format: file,m,1:3
sub check_fmt {
    my $in = shift(@_);
    my ($a, $b) = (1, 2);
    my $file = '';
    my @tabs = split(/\,/, $in);
    ($file, $a, $b) = @tabs;
    die("[$file] file not exists\n") if(! -e $file);
    die("[$in] error format: non-Integer found:\n") if(! $a =~ /^\d+$/ || ! $b =~ /^\d+|(\d+\:\d+)$/);
    if(@tabs > 3) {
        return 0;
    }else{
        if($b =~ /\:/) {
            die("[$in] not match FILE,N,m:m\n") if(! $b =~ /^\d+\:\d+$/);
            my ($n1, $n2) = $b =~ /^(\d+)\:(\d+)$/;
            die("[$in] not match FILE,N,m:m\n") if($n1 > $n2);
        }else {
            $b = $b . ':' . $b;
        }
        return ($file, $a, $b);
    }
}

sub usage {
    die(qq/
Usage: addNotes.pl [options] file1,1 file2,2 file3,2 ...

Add Notes from other files according to the id filed.

Options: -i     : The input txt file [id in column-1]
         -n     : The column contain the id field. [1]
         -o     : The output file

Note: the input file should be in the following format:

FILE,N,m

FILE: the input file
N   : 1 integer for the column of id if FILE
m   : INT, or INT:INT, to specify the note info in FILE

Example:
1. add info in col-7 of fileA to f.txt
perl addNotes.pl -i f.txt -o f.info fileA,1,7

2. add info in col-7 to col-10 of fileA to f.txt
perl addNotes.pl -i f.txt -o f.info fileA,1,7:10
\n/);
}

