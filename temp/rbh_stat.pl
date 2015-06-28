#!/usr/bin/env perl
use strict;
use warnings;

use File::Basename qw(basename);
use File::Spec::Functions qw(catfile);
use Getopt::Std;
use Data::Dumper;

#########################################
# Perform Reciprocal Best Hits analysis
#########################################

my $out = stat_mapping();
print $out . "\n";

##########################################
sub stat_mapping {
    my %opts = (f => 'ffn');
    getopts("d:o:", \%opts);
    usage() if(@ARGV != 1);
    die("[-d] Need input dir, contain *.ffn/*.faa\n") if(! defined $opts{d});
    die("[-o] Need input dir, contain blast output\n") if(! defined $opts{o});
    my $type = shift(@ARGV);
    die("[$type] unknown type: ffn or faa\n") if(! $type =~ /^(faa|ffn)$/);
    my $db_dir = $opts{d};
    my $hit_dir = $opts{o};
    #
    my @ids = glob("$db_dir/*.$type");
    @ids    = sort(@ids);
    my %stat = ();
    for(my $i = 0; $i < @ids; $i ++) {
        my $num1 = count_fasta($ids[$i]);
        my $name1 = basename($ids[$i]);
        $stat{$name1}->{$name1} = $num1;
        for( my $j = $i + 1; $j < @ids; $j ++) {
            my $name2 = basename($ids[$j]);
            my $best1 = catfile($hit_dir, $name1 . '_vs_'. $name2 . '.bestreciprocal');
            my $best2 = catfile($hit_dir, $name2 . '_vs_'. $name1 . '.bestreciprocal');
            my $num2 = count_tabline($best1, $best2);
            $stat{$name1}->{$name2} = $num2;
        }
    }
    my @gap = ();
    my $output = '';
    $output .= join("\t", 'ID', sort keys(%stat)). "\n";
    for my $a (sort keys %stat) {
        my @line = ();
        for my $b (sort keys %{$stat{$a}}) {
            push @line, $stat{$a}->{$b};
        }
        $output .= join("\t", $a, @gap, @line). "\n";
        push @gap, '-';
    }
    return $output;
}

#
sub count_fasta {
    my $in = shift(@_);
    my $num = 0;
    if( -e $in) {
        open my $fh_in, "< $in" or die "$!\n";
        while(<$fh_in>) {
            next unless(/^\>/);
            $num ++;
        }
        close $fh_in;
        return $num;
    }else {
        warn("[$in] file not found\n");
        return 0;
    }
}

sub count_tabline {
    my ($i, $j) = @_;
    my $num = 0;
    if( -e $i) {
        $num = read_file($i);
    }elsif( -e $j){    
        $num = read_file($j);
    }else {
        warn("[$i, $j] file not found\n");
    }
    return $num;
}

sub read_file {
    my $in = shift(@_);
    my $num = 0;
    if( -e $in) {
        open my $fh_in, "< $in" or die "$!\n";
        while(<$fh_in>) {
            next if(/^\s*$/);
            $num ++;
        }
        close $fh_in;
    }
    return $num;
}

sub usage {
    die("
Usage: stat_mapping.pl [options] <type>

Options: -d <STR>   : The dir of *.faa/*.fnn 
         -o <STR>   : The dir of blast output
         <type>     : [faa or ffn]

Example:
stat_mapping.pl -d db -o blast_out faa 
\n");
}


