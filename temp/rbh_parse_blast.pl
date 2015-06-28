#!/usr/bin/env perl

####################################################################
# Reciprocal Best Hits analysis using BLAST+
# 
# # bug, cannot process long_id sequence. (why, NC_002755/?)
####################################################################

use strict;
use warnings;
use File::Basename qw(basename dirname);
use File::Spec::Functions qw(catfile catdir);
use File::Path qw(make_path);
use Getopt::Std;
use Data::Dumper;

blast_rbh();
exit (1);

sub blast_rbh { # using blast+ for Reciprocal Best Hits analysis
    my %opts = (t => 'dna');
    getopts("o:n:a:b:h", \%opts);
    usage() if(defined $opts{h});
    die("[-o] Need specify the output dir\n") if(!defined($opts{o}));
    die("[-a, -b] Need input two fasta files\n") if(!defined $opts{a} || !defined($opts{b}));
    make_path($opts{o}) if(! -d $opts{o});

    my $best_outfile = catfile($opts{o}, $opts{n} . '.bestreciprocal');
    my $homo_file    = catfile($opts{o}, $opts{n} . '.homolist');
#    my $blast_ab = runBlast($opts{i}, $opts{d}, $opts{t}, $opts{o});
#    my $blast_ba = runBlast($opts{d}, $opts{i}, $opts{t}, $opts{o});
    my $blast_ab = processBlast($opts{a});
    my $blast_ba = processBlast($opts{b});    

    my $best_ab  = findBesthit($blast_ab);
    my $best_ba  = findBesthit($blast_ba);
    my ($best_pair, $homo_list) = bestReciprocal($best_ab, $best_ba);
    
    open my $fh_both, "> $best_outfile" or die "Cannot open $best_outfile, $!\n";
    open my $fh_homo, "> $homo_file" or die "Cannot open $homo_file, #!";
    print $fh_both join("\n", @$best_pair) . "\n";
    print $fh_homo join("\n", @$homo_list) . "\n";
    close $fh_both;
    close $fh_homo;
    print "Ouput RBH: " . scalar(@$best_pair) . "\n";
}

sub runBlast {
# BLAST+ 2.2.30+
# blastn/blastp: 
#   -outfmt 7 (tabular with comment lines)
    my ($query, $db, $type, $outdir) = (@_);
    die("[$query] not found\n") if(! -e $query);
    die("[$db] not found\n") if(! -e $db);
    my ($q_filename) = basename($query) =~ /(.*)\.\w+$/;
    my ($d_filename) = basename($db) =~ /(.*)\.\w+$/;
    my $blast_out = catfile($outdir, $q_filename . '_vs_' . $d_filename . '.out');
    # make blast db
    if( $type eq 'dna' ) {
        system "makeblastdb -dbtype nucl -parse_seqids -in $db -logfile $db.log";
        system "blastn -outfmt 7 -query $query -db $db > $blast_out ";
    }elsif( $type eq 'pro' ){
        system "makeblastdb -dbtype prot -parse_seqids -in $db -logfile $db.log";
        system "blastp -outfmt 7 -evalue 1e-2 -query $query -db $db > $blast_out"; # -use_sw_tback, very slow
    }
    my $blast_refine = processBlast($blast_out);
    return $blast_refine;
}

sub processBlast { # input blast+, -m7 output
    my $f = shift(@_); 
    open my $fh_f, "< $f" or die "Cannot open $f, $!\n";
    my $filt = '';
    while(<$fh_f>) {
        next if(m/^# [A-Z0]/); 
        $filt .= $_;
    }
    close $fh_f;
    return $filt;
}

# This version only select the first hit as best hit.
# waitting: add 1toN mode.
sub findBesthit {
    my $in = shift(@_);
    my @hits = split(/\#/, $in);
    shift(@hits);
    my %bt = ();
    for my $h (@hits) {
        my @hsp = split(/\n/, $h);
        shift(@hsp); # trim header 
# filt the blast output     
# @hsp = filtBlast(@hsp);
        my $rank = 0;
        for my $p (@hsp) {
            last if($rank < -3); ### only read 4 levels
            next if (! defined $p);
            my ($qid, $sid, $identity, $mismatch) = (split /\s+/, $p)[0, 1, 2, 4];
            next if($identity < 70 ); # discard low identity hits
#            next if($mismatch > 30); # number of mismatches
            next if(exists $bt{$qid}->{$sid}); # discard the following up duplicates
            $bt{$qid}->{$sid}->{'score'} = $rank;
            $bt{$qid}->{$sid}->{'info'}  = $p;
            $rank --;
        }
    }
    return (\%bt);
}

sub max_val {
    my $max = shift(@_);
    for my $n (@_){
        $max = ($max < $n)?$n:$max;
    }
    return $max;
}

sub min_val {
    my $min = shift(@_);
    for my $n (@_) {
        $min = ($min > $n)?$n:$min;
    }
    return $min;
}

sub bestReciprocal { # pass two hash: a->b, b->a, 
    my ($a, $b)  = (@_);
    my %hab      = %$a;
    my %hba      = %$b;
    my @best     = ();
    my @homolist = ();
# checking A->B
    for my $i (sort keys %hab) {
        my %check = ();
        for my $j (keys %{$hab{$i}}) {
            if( exists $hba{$j}->{$i} ) {
                $check{'sum'}->{$j} = $hab{$i}->{$j}->{'score'} +
                                      $hba{$j}->{$i}->{'score'};
                $check{'var'}->{$j} = abs($hab{$i}->{$j}->{'score'} -
                                          $hba{$j}->{$i}->{'score'});
            }
        }
# criteria
# max sum & min variance
        next if( keys %{$check{'sum'}} < 1 );
        my $max_sum = max_val( values %{$check{'sum'}} );
        my @hit_1s  = grep { $check{'sum'}->{$_} == $max_sum } keys %{$check{'sum'}}; # possible multiple hits with equal in score: 0, -2, and -1, -1
        my @vars    = ();
        for my $t (@hit_1s) {
            push @vars, $check{'var'}->{$t};
        }
        my $min_var = min_val( @vars );
        my @hit_2s  = grep { $check{'var'}->{$_} == $min_var } @hit_1s;
#
        if(@hit_2s > 1) {
            my $homo = join("\,", @hit_1s);
            push @homolist, '# possible homologs: ' . $i . "\t" . $homo;
            next;
        }
        push @best, $hab{$i}->{$hit_2s[0]}->{'info'};
        delete $hba{$hit_2s[0]};
        delete $hab{$i};
    }
# checking B->A
    for my $k (keys %hba) {
        my %check = ();
        for my $p (keys %{$hba{$k}}) {
            if( exists $hab{$p}->{$k} ) {
                $check{'sum'}->{$p} = $hba{$k}->{$p}->{'score'} +
                                      $hab{$p}->{$k}->{'score'};
                $check{'var'}->{$p} = abs($hba{$k}->{$p}->{'score'} -
                                          $hab{$p}->{$k}->{'score'});
            }
        }
        next if( keys %{$check{'sum'}} < 1);
        my $max_sum = max_val( values %{$check{'sum'}} );
        my @hit_1n  = grep { $check{'sum'}->{$_} == $max_sum } keys %{$check{'sum'}};
        my @vars    = ();
        for my $m (@hit_1n) {
            push @vars, $check{'var'}->{$m};
        }
        my $min_var = min_val( @vars );
        my @hit_2n  = grep { $check{'var'}->{$_} == $min_var } @hit_1n;

        if(@hit_2n > 1) {
            my $homo = join("\,", @hit_1n);
            push @homolist, '# possible homologs: ' . $k . "\t" . $homo;
            next;
        }
        push @best, $hab{$hit_2n[0]}->{$k}->{'info'};
    }
    return (\@best, \@homolist);
}

sub usage {
    die("
Usage: blast_rbh.pl [options] 

Options: -h     : show this info
         -o     : output dir
         -a     : output of blast+, with -m 7 format
         -b     : outout of blast+, with -m 7 format
         -n     : name prefix [default: out]

Examples:
1. find reciprocol best hit pairs between two fasta file

blast_rbh.pl -o out -n Ra2Rv -a a2b.blast -b b2a.blast 
\n");
}

##########
# change log
# 1. discard low identity hits (70%)
# 2. keep only the best pair of two ids.
# 3. evaluation: penalty by the ranking in each hit report. and variance.
#
#
#
#
#




