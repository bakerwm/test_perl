#!/usr/bin/env perl

#################################################
# Extract homologe sequences from a db 
# (multiple sequences)
#
# only output the first hit from blast output
# for each subject seq.
#
# 1. choose which name to output (NC_, name,...)
#
#################################################

use strict;
use warnings;
use File::Basename qw(basename dirname);
use File::Spec::Functions qw(catfile);
use File::Path qw(make_path remove_tree);
use Getopt::Std;

extract_homoseq_blast();
exit(1);

sub extract_homoseq_blast {
    my %opts = (n => 'id');
    getopts("d:n:o:f:", \%opts);
    usage() if(@ARGV != 1);
    die("[-d] Need input db file\n") if(! defined $opts{d});
    die("[-n $opts{n}] unknown -n\n") if(! $opts{n} =~ /^(id)|(name)$/);
    die("[-o] Need input output dir\n") if(! defined $opts{o});
    make_path($opts{o}) if(! -d $opts{o});
    my $query = shift(@ARGV);
    my $df = read_fa($opts{d});
    my $di = check_db_id($opts{d});
    my %db_fa = %$df;
    my %db_id = %$di;
# NC_123456, NC_123456.1, gi|123456|ref|NC_123456.1|, abcd efg chromosome, gi|...chromosome
    my $out_id = 0;
    $out_id    = 3 if($opts{n} eq 'name');
# run blast
    my $best_hit = run_blast($query, $opts{d});
    my %bt = %$best_hit;
    my $hit_info = '';
    for my $b (sort keys %bt) {
        my $q_file = catfile($opts{o}, $b . '.fa');
        open my $fh_q, "> $q_file" or die "Cannot open $q_file, $!\n";
        for my $t (sort keys %{$bt{$b}}) {
            die("[$t] not found in: $opts{d}\n") if(! exists $db_fa{$t});
            die("[$t] not found id: $opts{d}\n") if(! exists $db_id{$t});
            my ($start, $end) = (split /\t/, $bt{$b}->{$t})[8, 9];
            my $t_fa = pos2fa($start, $end, $db_fa{$t});
            my $t_id = $db_id{$t}[$out_id];
            print $fh_q join("\n", ">".$t_id, $t_fa)."\n";
            my $strand = '+';
            my ($bg, $ed) = ($start < $end)?($start, $end):($end, $start);
            my $length = $ed - $bg + 1;
            $hit_info .= $bt{$b}->{$t} . "\n";           
#$hit_info .= join("\t", $b, $t_id, $length, $bg, $ed, $strand). "\n";
        }
        close $fh_q;
    }
    # whether output info to file
    if(defined $opts{f} ) {
        open my $fh_f, "> $opts{f}" or die "Cannot open $opts{f}, $!\n";
        print $fh_f $hit_info;
        close $fh_f;
    }
    print "Finish: ouput [" . scalar(keys %bt) . '] files in [' . $opts{o} . "].\n";
}

sub check_db_id {
    my $in    = shift(@_);
    my %db_id = ();
    open my $fh_in, "< $in" or die "Cannot open $in, $!\n";
    while (<$fh_in>) {
        chomp;
        next if(! /^\>/);
# >gi|121635883|ref|NC_008769.1| Mycobacterium bovis BCG str. Pasteur 1173P2 chromosome, complete genome
        die("[$in] id is not standard NCBI format: $_\n") unless(/\>gi\|\d+\|\w+\|\w+\.\d+\|\s\w+/);
        $_ =~ s/\>//;
        my ($gi, $name) = split /\s/, $_, 2;
        my $nc_id    = (split /\|/, $gi)[3];
        my ($nc_id2) = $nc_id =~ /(\w+)\.\d+/;
        $name = s/\s/\_/;
        @{$db_id{$gi}} = ($nc_id2, $nc_id, $gi, $name, $_);
    }
    return(\%db_id);
}

#################################################################
# Perform blast search, return best hits
#################################################################
sub run_blast {
# BLAST+ 2.2.30+
# blastn/blastp:
# -outfmt 7 (tabular with comment lines)
    my ($query, $db) = (@_);
    die("[$query] not found\n") if(! -e $query);
    die("[$db] not found\n") if(! -e $db);
    my $q_filename = basename($query);
    my $d_filename = basename($db);
    my $rand_num = sprintf "%06d", int(rand(1000000));
    my $temp_dir = 'temp_' . $$ . '_' .  $rand_num;
    make_path($temp_dir);
    my $blast_out = catfile($temp_dir, $q_filename . '_vs_' . $d_filename . '.out');
    system "makeblastdb -dbtype nucl -in $db -logfile $db.log"; # do not use "-parse_seqids" for *.ffn
    system "blastn -outfmt 7 -query $query -db $db > $blast_out";
    my $hits = parsing_blast($blast_out);
    remove_tree($temp_dir);
    return $hits;
}

sub parsing_blast {
# parsing the blast output: -outfmt 7, tabular with comment lines
# report the best hit for each subject seq
    my $f = shift(@_);
    open my $fh_f, "< $f" or die "Cannot open $f, $!\n";
    my $filt = '';
    while(<$fh_f>) {
        next if(m/^# [A-Z0]/);
        $filt .= $_;
    }
    close $fh_f;
    #
    my @hits = split(/\#/, $filt);
    shift(@hits); # trim blank line
    my %besthit = ();
    for my $h (@hits) {
        my @hsp = split(/\n/, $h);
        shift(@hsp); # trim comment line
        for my $p (@hsp) {
            next if (! defined $p); # skip blank line
            my ($qid, $sid, $identity, $mismatch) = (split /\s+/, $p)[0, 1, 2, 4];
            next if($identity < 70 ); # discard low identity hits
#            next if($mismatch > 50 ); # too much mismatches
            next if(exists $besthit{$qid}->{$sid}); # only keep the first seq for each subject seq
            $besthit{$qid}->{$sid} = $p;           
        }
    }
    return(\%besthit);
}

################################
# Readin or output fasta seq
################################
sub read_fa {
    my $in   = shift(@_);
    my $text = '';
    open my $fh_in, "< $in" or die "Cannot open $in, $!\n";
    while(<$fh_in>) {
        $text .= $_;
    }
    close $fh_in;
    my @lines = split(/\>/, $text);
    shift(@lines); # trim blank line
    my %fa = ();
    for my $l (@lines) {
        my @g       = split(/\n/, $l);
        my $id_line = shift(@g);
        my $g_id    = (split /\s+/, $id_line)[0];
        my $g_fa    = join('', @g);
        $fa{$g_id} = $g_fa;
    }
    return(\%fa);
}

sub pos2fa {
    my ($start, $end, $line) = @_; 
    my $strand  = ($start < $end)?'+':'-';
    my ($s, $e) = ($start < $end)?($start, $end):($end, $start);
    my $length  = $e - $s + 1;
    my $begin   = $s - 1;
    my $fa = substr($line, $begin, $length);
    if($start > $end) {
        $fa =~ tr/ATCGatcg/TAGCtagc/;
        $fa = reverse($fa);
    }
    return($fa);
}

sub usage {
    die("
Usage: ext_homoseq.pl [options] <input.fa>

Options: -o <STR>   : Dir for output fasta seq
         -d <STR>   : Input a database file, in FASTA format
         -n <STR>   : How to name the new sequences?
                      id=NC_123456, name=M_tub_H37Rv [default: id]
         -f <STR>   : (optional), filename, write the position in in db

Example:
ext_homoseq.pl -o subseq -d db.fa -n id input.fa
\n");
}

