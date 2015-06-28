#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Std;

my %opts = ();
getopt("n:i:o:", \%opts);
if(!defined $opts{n} || !defined $opts{i} || !defined $opts{o}){
    &help();
    exit(1);
}
my @In_files = split(/\,/, $opts{i});
if($opts{n} != @In_files){
    print "The number of -n is not match the number of files in -i\n";
    exit(1);
}

sub help{
	print STDERR <<EOF;
Usage: perl merge4all.pl  -n 4  -i in_1,in_2,in_3,in_4  -o outfile
Options:
    -n <Integer>    : The number of input files.
    -i <file>	    : Input files.  
                      <1> File name in "45SE_H37Rv_*" format.                        
                      <2> Multiple files were delimited by comma ",".
    -o <file>	    : Output file.

EOF
}

#######################################################################
# Note: For multiple input files,
# 1. Filename like this, "45SE_H37Rv_*". eg: 45SE_H37Rv.txt
# 2. Multiple input filenames join by comma ",". eg: "infile_1,infile_2"
#
# Description:
# 1. Two sequences have an overlap longer than 40% of either of the 
#    sequence or longer than 10 nt will merge into a new one.
#    (match_length > 40% of input or 10 nt)
# 2. The expression level between sequences should not more than 
#    100-fold.
#    (Max_exp/exp <= 100)
#######################################################################

my $i = 0;
my @tag = ("a", "b", "c", "d");
my %hash;

foreach my $f (@In_files){
    open IN,$f or die;
    my $label = "n";
    if($f =~ /^\d+[S,P]E/i){
        $f =~ /^(\d+[S,P]E)\_/i;
        $label = $1;
    }
    while(<IN>){
        next if(/description/i);
        chomp;
        my @lines = split(/\t/,$_);
        $lines[0] = $label."-".$lines[0];
        @{$hash{$tag[$i]}->{$lines[0]}} = @lines;
    }
    close IN;
$i ++;
}

my $count_id = 0;
open O1,">$opts{o}" or die;
open O2,">$opts{o}\.match" or die;
foreach my $tag_n (sort keys %hash){
    foreach my $line_n (sort keys %{$hash{$tag_n}}){
        next unless(exists $hash{$tag_n}->{$line_n});
        my @subs = @{$hash{$tag_n}->{$line_n}};
        my ($ID, $EXP, $BE, $EN, $STR) = ($subs[0], $subs[1], $subs[3], $subs[4], $subs[5]);   
        $EXP = ($EXP =~ m/^\d+$|^\d+\.\d+$/)?$EXP:1;
        my ($min_begin, $max_end) = ($BE, $EN);
        my ($min_exp, $max_exp)   = ($EXP, $EXP);
        my %hash_out;
        foreach my $tag_m (sort keys %hash){
            foreach my $line_m (sort keys %{$hash{$tag_m}}){
                next unless(exists $hash{$tag_m}->{$line_m});
                my @querys = @{$hash{$tag_m}->{$line_m}};
                my ($id, $exp, $be, $en, $str) = ($querys[0], $querys[1], $querys[3], $querys[4], $querys[5]); 
                $exp = ($exp =~ m/^\d+$|^\d+\.\d+$/)?$exp:1;
# Same strand.
                next unless($STR eq $str);
                my $len   = $en - $be + 1;
                my $LEN   = $max_end - $min_begin + 1;
                my ($gap, $match) = (1,1);
                if($be < $min_begin && $en > $min_begin){
                    $gap   = $en - $min_begin; 
                    $match = ($gap>$LEN)?$LEN:$gap;
                }elsif($be >= $min_begin && $be < $max_end){
                    $gap = $EN - $be; 
                    $match = ($gap>$len)?$len:$gap;
                }else{
                    next;
                }
# overlap is 40% of either of the sequences or larger than 10 bp.
                next unless($match/$len >= 0.4 || $match/$LEN >= 0.4 || $match >= 10);
                my $query      = join("\t",@querys); 
                $hash_out{$id} = $query;
                ($min_exp, $max_exp)   = &min_max($min_exp, $max_exp, $exp);
                ($min_begin, $max_end) = &min_max($min_begin, $max_end, $be, $en);
                delete $hash{$tag_m}->{$line_m};
            }
        }
        $count_id ++;
        my $print_out = "";
        my ($final_begin, $final_end) = (10000000, 0);
        foreach(keys %hash_out){
            my @temp = split (/\t/, $hash_out{$_});
            my ($t_exp, $t_begin, $t_end) = ($temp[1], $temp[3], $temp[4]);
            $t_exp = ($t_exp =~ m/^\d+$|^\d+\.\d+$/)?$t_exp:1;
# No more than 100-fold exp between merged sequences.
            next if($max_exp/$t_exp > 100);
            $print_out .= "$hash_out{$_}\n";
            ($final_begin, $final_end) = &min_max($final_begin, $final_end, $t_begin, $t_end);
        }
        my $out_id    = sprintf "%.4d", $count_id;
        my $final_len = $final_end - $final_begin + 1;
        my $seed_out  = join"\t",("Seed$out_id", $max_exp, $final_len, $final_begin, $final_end, $STR);
        my $match_out = join"\t",(">Seed$out_id",$max_exp,  "$final_begin\:$final_end\:$STR");
        print O1 $seed_out,"\n";
        print O2 $match_out,"\n",$print_out;
        undef %hash_out;
    }
}
close O1;
close O2;

sub min_max{
    my $min = shift(@_);
    my $max = shift(@_);
    foreach (@_){
        $min = ($min < $_)?$min:$_;
        $max = ($max > $_)?$max:$_;
    }
    return ($min, $max);
}
