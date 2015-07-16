#!/usr/bin/env perl

#########################################################
# It is a converter that will convert the txt files to 
# plot/histogram data/labels for Circos
#
# # Plot data (4-column)
# chr1 Start End Value 
# 
# # label data
# In Value field.
#########################################################

use strict;
use warnings;
use File::Basename qw(basename dirname);
use File::Path qw(make_path remove_tree);
use Getopt::Std;

usage() if(@ARGV == 0);
my $command = shift(@ARGV);
my %prog    = (plot => \&sort2circosplot,
               text => \&sort2label,
               link => \&sort2link);
die("Unknown command [$command] \n") if(! defined($prog{$command}));

&{$prog{$command}};

exit(0);

# sort to plots data for circos
sub sort2circosplot {
    my %opts = (s => '+/-', k => 1, m => 1);
    getopts('o:s:k:a:m:c', \%opts);
    die(qq/
Usage: sort2circos.pl plot [options] <input.txt>

Options: -o STR     The dir for the output
         -c STR     Add color to the 4-column of output
         -s STR     The strand of file: +, - or +\/-, [+\/-]
         -k INT     Value in k-column to present in output file.[1]
         -a STR     The annotation data for the input (expression, color, ...)
                    col-1 should match the col-1 of input file.
         -m INT     Value in N-column of the annotation file [1]

Example: 
1. Present column-3 of infile to output:
sort2circos.pl plot -k 3 infile.txt

2. Add column-8 of annotation file to output:
sort2circos.pl plot -a anno.txt -m 8 infile.txt

3. Add random-color to value-column
sort2circos.pl plot -c infile.txt 

4. Convert strand + from infile
sort2circos.pl plot -s + infile.txt 
\n/) if(@ARGV == 0);

#    die("[-o] Need specify the output dir\n") if(! defined $opts{o});
    die("[-s $opts{s}] need to be +, -, +/-") unless($opts{s} =~ /^[+-]$|^\+\/\-$/);
    die("[-k $opts{k}] id expected a number\n") if(! $opts{k} =~ /^\d+$/);
    die("[-m $opts{m}] is expected a number\n") if(! $opts{m} =~ /^\d+$/);
    my $infile = shift(@ARGV);
    die("[$infile] is empty\n") if(! not_empty_file($infile) );

    # read txt
    my %info = ();
    my %anno = ();
    %info    = readinfo($infile,  $opts{k}); # return the whole line
    %anno    = readinfo($opts{a}, $opts{m}) if( defined $opts{a});

    # output format
    open my $fh_in, "< $infile" or die "Cannot open $infile, $!\n";
    while(<$fh_in>) {
        next if(/(^\s*$)|(^\#)/);
        chomp;
        my @tabs = split("\t", $_);
        my $val  = $info{$tabs[0]};
        if(exists $anno{$tabs[0]}) {
            $val = $anno{$tabs[0]};
        }
        # select strand
        my $strand = $tabs[5];
        next if(index($opts{s}, $strand) == -1);
        # select color
        if (defined $opts{c}) {
            $val = randcolor();
        }
        print join("\t", 'chr1', $tabs[3], $tabs[4], $val) . "\n";
    }
    close $fh_in;
}

#
sub not_empty_file {
    my $in = shift(@_);
    my $count = 0;
    if(-e $in) {
        open my $fh_in, $in or die "Cannot open $in, $!\n";
        while(<$fh_in>) {
            chomp;
            next if(/(^\s*$)|(^\#)/);
            $count ++
        }
        close $fh_in;
    }
    return $count;
}

sub readinfo {
    my ($in, $col)   = @_;
    my %info = ();
    if( not_empty_file($in) ) {
        open my $fh_in, $in or die "Cannot open $in, $!\n";
        while(<$fh_in>) {
            chomp;
            next if(/(^\s*$)|(^\#)/);
            my @ts = split("\t", $_);
            die("[-k|-m] -column not found\n") if(@ts < $col);
            $info{$ts[0]} = $ts[$col - 1];
            if($col == 0) {
                $info{$ts[0]} = '';
            }
        }
        close $fh_in;
    }
    return %info;
}


sub randcolor {
    my $i = int(rand(22));
    return 'fill_color=chr'.$i;
}

sub usage {
    die(qq/
Usage: sort2circos.pl <command> [options] \n
Command: plot   Change the input txt to circos plots: histogram, heatmap, 
                line, tile, ...
         text   Change the input txt to text labels:
         line   Change the input txt to link plots:
\n/);
}

