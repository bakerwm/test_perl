#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Std;

#########################################
# Perform Reciprocal Best Hits analysis
#
# input db/
#########################################

my $blast_rhb = '/home/wangming/work/bin/temp/rbh_run_blast.pl';

run_pairs();
exit(1);

## subroutines
sub run_pairs {
    my %opts = (p => 4);
    getopts("i:o:p:", \%opts);
    usage() if(@ARGV != 1);
    die("[-i] input the db dir:\n") if(! defined $opts{i});
    die("[-o] output of blast files: \n") if(! defined $opts{o});
    die("[-p] need to be INT [2-20]\n") if(! $opts{p} =~ /^\d+$/);
    my $type = shift(@ARGV);
    die("[$type] unknown type: ffn, faa\n") if(! $type =~ /^(ffn|faa)$/);
    my @ids = glob("$opts{i}/*.$type");
    @ids = sort(@ids);
    my @runs = ();
    for(my $i = 0; $i < @ids; $i ++) {
        for(my $j = $i + 1; $j < @ids; $j ++) {
            my $cmd = run_rbh($ids[$i], $ids[$j], $type, $opts{o});
            push @runs, $cmd;
        }
    }
    mp_sys($opts{p}, \@runs);
}

sub run_rbh {
    my ($a, $b, $t, $out) = @_;
    my $name = ($t eq 'ffn')?'dna':(($t eq 'faa')?'pro':'-');
    my $cmd = "$blast_rhb -o $out -t $name -i $a -d $b ";
    return $cmd;
}

sub mp_sys {
    my ($max, $c) = @_;
    my @commands = @$c;
    my $children = 0;
    my @pids = ();
    for(my $i = 0; $i < @commands; $i ++) {
        my $pid;
        if($children >= $max) {
            $pid = wait();
            $children --;
        }
        $pid = fork();
        die("Cannot fork\n") if (! defined $pid);
        if( $pid ){
            $children ++;
            push @pids, $pid;
        }else {
            system("$commands[$i]");
            exit 0;
        }
    }
    for my $n (@pids) {
        my $chk = waitpid($n, 0);
    }
    sleep(5);
}

sub usage {
    die("
Usage: run_pairs.pl [options] <type>

Options: -i <STR>   : Dir of fasta files, *.faa, *.ffn
         -o <STR>   : Dir of blast output files
         -p <INT>   : Number of children processes at once, [4]
         <type>     : faa or ffn

Example:
run_pairs.pl -i db -o output -p 10 ffn
\n")
}

