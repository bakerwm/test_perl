#!/usr/bin/env perl

############################################################
# Search the Rfam database for possible sRNAs (input.fa)
#
# Criteria:
# 1. E-value < 0.01 (col 17 is '!')
#
# # Step 1.
# cmscan with default parameters:
# cmscan --tblout hit.tbl Rfam.cm input.fa > input_Rfam.out
# cat hit.tbl | grep -v \# | grep \! > hit.filterout.tbl
# 
# # Step 2.
# parse the tbl output
# query_id  hit_id1,hit_id2,...
##############################################################

use strict;
use warnings;
use File::Basename qw(basename dirname);
use File::Spec::Functions qw(catfile catdir);
use File::Path qw(make_path remove_tree);
use File::Which qw(which);
use POSIX qw(strftime);
use Getopt::Std;

searchRfam();
exit (1);

###
sub searchRfam {
    my %opts = (d => '/share/wangming/database/Rfam/current/Rfam.cm');
    getopts("o:f:d:s:", \%opts);
    usage() if (@ARGV == 0);
    my $input = shift(@ARGV);
    die("<$input> file not exists\n") if(! -e $input);
    my $cmscan = 'cmscan';
    if(defined $opts{s}) {
        $cmscan = $opts{s};
        die("[-s $opts{s}] file not exist\n") if(! -e $opts{s});
    }else{
        die("[cmscan] not found. use -s to specify the path of cmscan\n") 
            if(! which('cmscan'));
    }
    # specific the output dir
    die("[-o $opts{o}] need specify the output dir\n") if(!defined $opts{o});
    make_path($opts{o}) if(! -d $opts{o});

    # perform cmscan
    print show_date('Running cmscan') . "\n";
#    print '[' . show_date() . ']' . "\t" . 'Running cmscan' . "\n"; 
    my ($in_prefix) = basename($input) =~ /(.*)\.\w+$/;
    my $Rfam_out    = catfile($opts{o}, $in_prefix . '_Rfam.out');
    my $Rfam_tbl    = catfile($opts{o}, $in_prefix . '_Rfam.tbl');
    my $Rfam_filt   = catfile($opts{o}, $in_prefix . '_Rfam.filt.tbl');
    my $Rfam_index  = catfile($opts{o}, $in_prefix . '_Rfam.index');
    system "$cmscan --tblout $Rfam_tbl $opts{d} $input > $Rfam_out";
#    my $run_search  = "$cmscan --tblout $Rfam_tbl $opts{d} $input > $Rfam_out";
#    my $run_file    = 'cat '. $Rfam_out . ' | grep -v \# | grep \! > ' . $Rfam_filt;

    # wrap output
    # criteria : !
    # output format: acc:id,acc:id,...
    print show_date('Wrap results') . "\n";
#    print '[' . show_date() . ']' . "\t" . 'Wrap results' . "\n";
    open my $fh_in, "< $Rfam_tbl" or die "Cannot open $Rfam_tbl, $!\n";
    open my $fh_out, "> $Rfam_index" or die "Cannot open $Rfam_index, \n";
    my %hit = ();
    while(<$fh_in>) {
        chomp;
        next if(/\#|(^\s*$)/); # skip blank lines
        next if(! /\!/); # indicates whether or not this hit achieves the inclusiv threshold. '!'
        my @tabs = split(/\s+/, $_);
        push @{$hit{$tabs[2]}}, $tabs[1] . ':' . $tabs[0];
    }
    close $fh_in;
    my %anno = ();
    for my $h (sort keys %hit) {
        $anno{$h} = join("\,", @{$hit{$h}});
        print $fh_out join("\t", $h, $anno{$h}) . "\n";
    }
    close $fh_out;

    # add annotation to txt file
    my $txt_anno   = catfile($opts{o}, $in_prefix . '.anno.txt');
    if(defined($opts{f})) {
        die("[-f $opts{f}] file not exists\n") if(! -e $opts{f});
        annoTxt(\%anno, $opts{f}, $txt_anno);    
    }
    print show_date('Finish') . "\n";
    print 'Find results in [' . $opts{o} . ']' . "\n";
}

sub annoTxt {
    my ($info, $filein, $fileout) = @_;
    my %hinfo = %{$info};
    open my $fh_in, "< $filein" or die "Cannot open $filein, $!\n";
    open my $fh_out, "> $fileout" or die "Cannot open $fileout, $!\n";
    while(<$fh_in>) {
        chomp;
        if(/\#|(^\s*$)/) {
            print $fh_out $_ . "\n";
        }
        my @tabs = split(/\s+/, $_);
        my $note = (exists $hinfo{$tabs[0]})?$hinfo{$tabs[0]}:'-';
        print $fh_out join("\t", @tabs, $note) . "\n";
    }
    close $fh_in;
    close $fh_out;
}

sub show_date {
    my $text = shift(@_);
    my $date = strftime "%Y-%m-%d %H:%M:%S", localtime;
    my $log  = '['. $date . ']' . "\t" . $text;
    return $log;
}

sub usage {
    die(qq/
Usage: searchRfam.pl [options] <in.fa>

Options: -o     : The dir for output files
         -f     : the original tab file of input
         -d     : path to the Rfam.cm database
                  [default: \/share\/wangming\/database\/Rfam\/current\/Rfam.cm]
         -s     : path to command: [default: cmscan]

Examples:
1. only output Rfam output (id)
searchRfam.pl -o out input.fa

2. add Rfam annotation to a txt file
searchRfam.pl -o out -f input.txt -d Rfam.cm input.fa
\n/);
}
