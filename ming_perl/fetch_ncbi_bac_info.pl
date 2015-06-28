#!/usr/bin/env perl

##################################################
# Download the LIST of bacteria genomes in NCBI
# from the FTP site:
# ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria
#
# Usage: perl ncbi_bac_genome_list.pl  
#
# use `GET` command to fetch content of a url
#
# WangMing wangmcas(AT)gmail.com
# 2015-06-22 
##################################################

use strict;
use warnings;
use Getopt::Std;
use POSIX qw(strftime);

my %genomes = ( 1 => 'Bacteria',
                2 => 'Fungi');

fetch_ncbi_acc();
exit(1);

# write the bacteria genome names to file: 
#
# Example:
# ncbi_bacteria_genomes_20150622_14-20-00.txt
#
sub fetch_ncbi_acc {
    my %opts = (p => 10, n => 1);
    getopts("p:n:h", \%opts);
    die("[-p|-n] Input Integer:\n") unless($opts{p} =~ /^\d+$/ && $opts{n} =~ /^\d+$/);
    die("[-n] unknown command: 1 or 2\n") unless(defined($genomes{$opts{n}}));
    usage() if(@ARGV != 1);
    my $date = strftime "%Y%m%d", localtime;
    my $outfile = $ARGV[0] . '_' . $date . '.txt';
    my $max = $opts{p};
    my $url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/' . $genomes{$opts{n}};
    print '['. show_date() . ']'. ' Parsing genome names:' . "\n";
    # main process
    my @ids = parse_genome_id($url);
    print '['. show_date() . ']'. ' Parsing accession IDs for each genome:'. "\n";
    my @out = mp_runs(\@ids, $outfile, $max, $url); # using fork to clone N process
}

sub show_date {
    my $date = strftime "%Y-%m-%d %H:%M:%S", localtime;
    return $date;
}

# start N child processes at once ()
#
sub mp_runs {
    my ($n, $f, $max, $url) = @_;
    my @names = @$n;
    my @pids  = ();
    my $children = 0;
    open my $fh_f, "> $f" or die "Cannot open $f, $!\n";
    for(my $i = 0; $i < @names; $i++) {
        my $position = $i + 1;
        print '['. show_date() .'] - ' . $names[$i] . ' ['. $position .'/'. scalar(@names) .']' . "\n";
        my $pid;
        if($children >= $max) {
            $pid = wait();
            $children --;
        }
        # clone a child process
        $pid = fork(); # clone a proc (child)
        die("Cannot fork\n") if(! defined $pid);
        if($pid) { # parent proc
            $children ++; # add a child process
            push @pids, $pid; # 
        }else { # this is child proc
            my $info = child($names[$i], $url); # run in child proc
            print $fh_f $info . "\n";
            exit 0; # terminate this child proc
        }
    }
    # wait for all child to finish.
    for my $n (@pids) {
        my $chk  = waitpid($n, 0);
        my $info = $? >> 8; # remove signal / dump bits from rc
        print "PID $n finished with info $info\n";
    }
    close $fh_f;
    # run child proc
    sub child {
        my ($id, $url) = @_;
        my $id_acc  = read_genome_dir($id, $url);
        sleep(2);
        return "$id\t$id_acc";
    }
}

# parsing the genome IDs
#
# Examples:
# dr-xr-xr-x   2 ftp      anonymous     4096 Dec  6  2010 Mycobacterium_tuberculosis_H37Rv_uid57777
#
sub parse_genome_id {
    my $url = shift(@_);
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
    my ($id, $url) = @_;
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

# usage
sub usage {
    die(qq/
Usage: fetch_ncbi_info.pl [Options] <output>

Options: -p <INT>   : start N child process at once: (2-20) [10]
         -n <INT>   : choose genomes: 1=Bacteria, 2=Fungi [1]
         output     : write the result to file

Example:
fetch_ncbi_info.pl -p 20 -n 1 Bacteria > bac.log
fetch_ncbi_info.pl -p 10 -n 2 Fungi    > fungi.log
\n/);
}
