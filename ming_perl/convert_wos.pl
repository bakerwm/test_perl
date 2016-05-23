#!/usr/bin/perl -w

##-------------------------------------------------------------------------
## description
##
## This script was designed to process the output of Web of Science records
## and convert the output to tab-separated format.
## please refer http://webofknowledge.com for more details.
##
## Get the input WOS file:
## "Web of Science Core Collection" -> "Hit search" -> "Save to Other 
## File Formats" -> "Record Content: Full Record and Cited References" & 
## "File Format: Plain text"
##
## Ming Wang (wangmcas@gmail.com)
## 2016-05-23
##-------------------------------------------------------------------------

use strict;
use warnings;
use Getopt::Std;
use POSIX qw(strftime);

sub show_time {
    my $date = strftime "%Y-%m-%d %H:%M:%S", localtime;
    print STDERR "[$date] - ". $_[0] . "\n";
}

sub usage {
    die"
Usage: perl $0 [-f|-d] <WOS.txt>

Options:
    -f <STR> :		Selected fields of formats, connected by comma,
                        defualt: [PY,AU,TI,J9,DI,EM]
    -d       :		if specified, only report the duplicate records
\n";
exit;
}

##---------------------------
## parsing parameters
##---------------------------
my %opts = ();
getopts("f:d", \%opts);
&usage if(@ARGV == 0 && -t STDIN);
my $is_dup = defined($opts{d});
my $in_fmt = 'PY,AU,TI,J9,DI,EM';
$in_fmt = $opts{f} if defined($opts{f});
my @in_fmt_tags = split ",", $in_fmt;

##---------------------------
## parsing input file
##---------------------------
my $counter = 0;
my $flag = "NULL";
my $var  = "NULL";
my %ha = ();
my %flags = ();
&show_time("Start parsing WOS records");
while(<>) {
    chomp;
    next if(/^\s+$/);
    $counter ++ if(/^PT\s/);
    if(/^\w+/) {
        ($flag, $var) = split "\\s", $_, 2;
	$var = (defined $var)?$var:"NULL";
	$flags{$flag} ++;
    }else {
        $var = $_;
	$var =~ s/^\s+//g;
    }
    push @{$ha{$counter}->{$flag}}, $var;
}

##---------------------------
## check supported formats
##---------------------------
my $fmt_status = 1;
my @fmt_checks = ();
for my $f (@in_fmt_tags) {
    if(! exists $flags{$f}) {
        $f .= "*";
	$fmt_status *= 0;
    }
    push @fmt_checks, $f;
}
if($fmt_status == 0) {
    print STDERR join(",", @fmt_checks) . "\n" . '[* field(s) not recognized]' . "\n";
    exit;
}

##---------------------------
## convert formats
##---------------------------
my @raw_output = ();
my %uniq = ();
for my $n (sort keys %ha) {
    my @out_fmts = ();
    for my $fmt (@in_fmt_tags) {
        my $var = (defined $ha{$n}->{$fmt})?(join " ", @{$ha{$n}->{$fmt}}):"NULL";
	push @out_fmts, $var;
    }
    my $out_fmt = join("\t", @out_fmts);
    next if($out_fmt =~ /^NULL/);
#    my $year    = (defined $ha{$n}->{"PY"}[0])?$ha{$n}->{"PY"}[0]:"NULL";
#    my $author  = (defined $ha{$n}->{"AU"}[0])?$ha{$n}->{"AU"}[0]:"NULL";
#    my $title   = (defined $ha{$n}->{"TI"}[0])?(join " ", @{$ha{$n}->{"TI"}}):"NULL";
#    my $journal = (defined $ha{$n}->{"J9"}[0])?$ha{$n}->{"J9"}[0]:"NULL";
#    my $doi     = (defined $ha{$n}->{"DI"}[0])?$ha{$n}->{"DI"}[0]:"NULL";
#    my $email   = (defined $ha{$n}->{"EM"}[0])?$ha{$n}->{"EM"}[0]:"NULL";
#    my $out_fmt = join("\t", $year, $author, $title, $journal, $doi, $email);
    push @raw_output, $out_fmt;
    $uniq{$out_fmt} ++;
}

##---------------------------
## output
##---------------------------
my $num = 0;
for my $fmt (sort keys %uniq) {
    if($is_dup) {
        if($uniq{$fmt} > 1) {
	    print $fmt . "\n";
	    $num ++;
	}
    }else {
        print $fmt . "\n";
	$num ++;
    }
}

&show_time("Processed records: $num");
&show_time("Finish!");

## END ##
