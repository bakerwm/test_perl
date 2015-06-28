#!/usr/bin/env perl

##################################################
# Download the LIST of bacteria genomes in NCBI
# from the FTP site:
# ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria
#
# Usage: perl ncbi_bac_genome_list.pl  
# # waiting for XX minutes, depends on your network.
#
# WangMing wangmcas(AT)gmail.com
# 2015-06-22 
##################################################

use strict;
use warnings;
#use LWP::Simple;
use POSIX qw(strftime);

my $url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria';
print "Usage: perl ncbi_bac_genome_lists.pl\n";
fetch_ncbi_acc();
exit(1);

# write the bacteria genome names to file: 
#
# Example:
# ncbi_bacteria_genomes_20150622_14-20-00.txt
#
sub fetch_ncbi_acc {
    print '['. show_date() . ']'. ' Parsing genome names:' . "\n";
    my @ids = parse_genome_id();
    print '['. show_date() . ']'. ' Parsing accession IDs for each genome:'. "\n";
    my $date = strftime "%Y%m%d", localtime;
    my $bak_file = 'ncbi_bacteria_genomes_' . $date . '.txt';
    open my $fh_bak, "> $bak_file" or die "Cannot open file $bak_file, $!\n";
    my $count = 1;
    for my $g (sort @ids) {
        print '['. show_date() .'] - ' . $g . ' ['. $count .'/'. scalar(@ids) .']' . "\n";
        my $id_acc  = read_genome_dir($g);
        $count ++;
        print $fh_bak $g . "\t" . $id_acc . "\n";
    }
    close $fh_bak;
}

sub show_date {
    my $date = strftime "%Y-%m-%d %H:%M:%S", localtime;
    return $date;
}

sub mp_runs {
    my @commands = @_;
    my @pids;
    for my $cmd ( @commands ) {
        my $pid = fork;
        if( $pid ) {
            push @pids, $pid;
            next;
        }
#        system( $cmd );
        my $id_acc = read_genome_dir($g);
        exit;
    }
    wait for @pids;
    system "sleep 5";
}

# parsing the genome IDs
#
# Examples:
# dr-xr-xr-x   2 ftp      anonymous     4096 Dec  6  2010 Mycobacterium_tuberculosis_H37Rv_uid57777
#
sub parse_genome_id {
    my $content = qx{GET $url};
    die "Could not get $url\n" unless defined $content;
    my @lines = split /\n/, $content;
    # parsing the genome IDs
    my @genomes = ();
    for(@lines) {
        chomp;
        next unless(my ($name) = $_ =~ /\d+\s+(\w+)$/); # one strain in each subdirectory
        push @genomes, $name;
    }
    return @genomes;
}

# parsing the content of each genome
#
# Examples:
# -r--r--r--   1 ftp      anonymous   157421 Oct 22  2010 NC_009932.fna
#
sub read_genome_dir {
    my $id = shift(@_);
    my $sub_url = $url . '/' . $id;
    my $sub_content = qx{GET $sub_url};
    warn "Could not get $sub_url\n" unless defined $sub_content;
    my @sub_lines = split /\n/, $sub_content;
    my @acc_ids = ();
    for(@sub_lines) {
        chomp;
        next unless(my ($acc) = $_ =~ /\d+\s+(NC\_\d+)\.fna/); # parsing the *.fna files
        push @acc_ids, $acc;    
    }
    my $note = join("\,", @acc_ids);
    return $note;
}

