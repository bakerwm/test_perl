#!/usr/bin/env perl

######################################
# compare lines between a anf b
#
######################################

use strict;
use warnings;
use File::Basename qw(basename dirname);
use File::Spec::Functions qw(catfile catdir);
use Getopt::Std;
use Data::Dumper;

my %opts = ();
$opts{m} = $opts{n} = 1;
getopt("m:n:o:",\%opts);
die "Usage: a2b.pl [-m] <a-col> [-n] <b-col> [-o] <out-dir> fileA fileB\n" if(@ARGV != 2);

my $file_a = shift;
my $file_b = shift;

my ($stat, $a_only, $b_only, $a_both, $b_both) = 
    &cmp_both_files($file_a, $file_b, $opts{m}, $opts{n});

print $stat,"\n";

if($opts{o}) {
    my $outdir = $opts{o};
    mkdir $outdir unless -d $outdir;
    my ($a_pre, $a_ext) = basename($file_a) =~ /(.*)\.(\w+)$/;
    my ($b_pre, $b_ext) = basename($file_b) =~ /(.*)\.(\w+)$/;
    my $file_a_both = catfile($outdir, "$a_pre\_both\.$a_ext");
    my $file_a_only = catfile($outdir, "$a_pre\_only\.$a_ext");
    my $file_b_both = catfile($outdir, "$b_pre\_both\.$b_ext");
    my $file_b_only = catfile($outdir, "$b_pre\_only\.$b_ext");
    
    &write_file($a_both, $file_a_both);
    &write_file($a_only, $file_a_only);
    &write_file($b_both, $file_b_both);
    &write_file($b_only, $file_b_only);
}


# Subroutines #
sub write_file{
    my ($var, $file) = @_;
    open OUT, "> $file" or die "$!";
    print OUT $var,"\n";
    close OUT;
}

sub cmp_both_files{
    my ($a, $b, $col_a, $col_b) = @_;
    my %id_a = &readfile($a, $col_a);
    my %id_b = &readfile($b, $col_b);
    # check a to b 
    my @A_both = ();
    my @B_both = ();
    my @A_only = ();
    my @B_only = ();
    foreach my $ta (keys %id_a) {
        if(exists $id_b{$ta}) {
            push @A_both, $id_a{$ta};
        } else {
            push @A_only, $id_a{$ta};
        }
    }
    # check b to a
    foreach my $tb (keys %id_b) {
        if(exists $id_a{$tb}) {
            push @B_both, $id_b{$tb};
        } else {
            push @B_only, $id_b{$tb};
        }
    }
    # stat num
    my $len_a = keys %id_a;
    my $len_b = keys %id_b;
    my $len_both = @A_both;
    my $len_a_only = $len_a - $len_both;
    my $len_b_only = $len_b - $len_both;

    my $a_name = basename($a);
    my $b_name = basename($b);
    my $stat_out = join "\n", ("\#\t$a_name\t$b_name", 
                               "Total:\t$len_a\t$len_b", 
                               "Both:\t$len_both\t$len_both",
                               "Uniq:\t$len_a_only\t$len_b_only"
                               );
    my ($file_a_only, $file_a_both, $file_b_only, $file_b_both);
    $file_a_only = join "\n", (sort @A_only);
    $file_b_only = join "\n", (sort @B_only);
    $file_a_both = join "\n", (sort @A_both);
    $file_b_both = join "\n", (sort @B_both);
    return ($stat_out, $file_a_only, $file_b_only, $file_a_both, $file_b_both);
}

sub readfile{
    my ($in, $col) = @_;
    my %id = ();
    open F, "< $in" or die "$!";
    while(<F>){
        chomp;
        next if(/(^\s*$)|(^\#)/);
        my $col_index = $col - 1;
        my @tabs = split /\s+/;
        my $name = $tabs[$col_index];
        push @{$id{$name}}, $_;
    }
    close F;
    # check id replicates
    my %hn = ();
    foreach my $n (keys %id) {
        die "Error: found replicates in column [$col] of [$in]\n" if(@{$id{$n}} > 1);
        $hn{$n} = shift(@{$id{$n}});
    }
    return %hn;
}

