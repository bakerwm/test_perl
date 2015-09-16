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
    my %opts = (t => 'NC_000962', n => 1);
    getopts("d:t:n:o:f:", \%opts);
    usage() if(@ARGV != 1);
    die("[-d] Need input db file\n") if(! defined $opts{d});
    die("[-n $opts{n}] unknown -n\n") if(! $opts{n} =~ /^(id)|(name)$/);
    die("[-o] Need input output dir\n") if(! defined $opts{o});
    make_path($opts{o}) if(! -d $opts{o});
    my $query = shift(@ARGV);
    my %db_fa = %{ read_fa($opts{d}) };
    my %db_id = %{ check_db_id($opts{d}) };
# 0=NC_123456, 1=NC_123456.1, 2=gi|12344678|ref|NC_123456.1|, 3=gi|....chromosome
    my %id_type = (1 => 0, 2 => 2, 3 => 3 );
    if(! exists $id_type{$opts{n}}) {
        die("[$opts{n}] unknown type, 1, 2 or 3\n");
    }
    my $id_fmt = $id_type{$opts{n}};
# run blast
    my %best_hit = %{ run_blast($query, $opts{d}) };
    my $hit_info = '';
    for my $q (sort keys %best_hit) { # query
        my %q_fa = ();
        for my $s (sort keys %{$best_hit{$q}}) {
            die("[$s] not found in: $opts{d}\n") if(! exists $db_fa{$s});
            die("[$s] not found id: $opts{d}\n") if(! exists $db_id{$s});
            my ($start, $end) = (split /\t/, $best_hit{$q}->{$s})[8, 9];
            my $s_fa = pos2fa($start, $end, $db_fa{$s});
            my $s_id = $db_id{$s}[$id_fmt];
            $q_fa{$s_id} = $s_fa;
            my $strand = '+';
            my ($bg, $ed) = ($start < $end)?($start, $end):($end, $start);
            my $length = $ed - $bg + 1;
            $hit_info .= $best_hit{$q}->{$s} . "\n";
#$hit_info .= join("\t", $b, $t_id, $length, $bg, $ed, $strand). "\n";
        }
        ### chagnge the output format
        my $q_file = catfile($opts{o}, $q . '.fa');
        open my $fh_q, "> $q_file" or die "Cannot open $q_file, $!\n";
        if(exists $q_fa{$opts{t}}) {
            print $fh_q join("\n", '>'.$opts{t}, $q_fa{$opts{t}} ). "\n";
            delete $q_fa{$opts{t}};
            for my $f (keys %q_fa) {
                print $fh_q join("\n", '>'.$f, $q_fa{$f}). "\n";
            }
        }else{
            die "[$q] target ref: $opts{t} not found, check (-n, -t)\n";
        }
        close $fh_q;
    }
    # whether output info to file
    if(defined $opts{f} ) {
        open my $fh_f, "> $opts{f}" or die "Cannot open $opts{f}, $!\n";
        print $fh_f $hit_info;
        close $fh_f;
    }
    print "Finish: ouput [" . scalar(keys %best_hit) . '] files in [' . $opts{o} . "].\n";
}

sub check_db_id {
    my $in    = $_[0];
    my %db_id = ();
    open my $fh_in, "< $in" or die "Cannot open $in, $!\n";
    while (<$fh_in>) {
        chomp;
        next if(! /^\>/);
# >gi|121635883|ref|NC_008769.1| Mycobacterium bovis BCG str. Pasteur 1173P2 chromosome, complete genome
        die("[$in] id is not standard NCBI format: $_\n") unless(/\>gi\|\d+\|\w+\|\w+\.\d+\|\s\w+/);
        $_ =~ s/\>//;
        my ($gi, $name) = split /\s/, $_, 2;     
        my $nc_id       = (split /\|/, $gi)[3];  # NC_000962.3
        my ($nc_id2)    = $nc_id =~ /(\w+)\.\d+/; # NC_000962
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
         -n <INT>   : How to name the new sequences?
                      1=NC_123456, 2='gi|123456|ref|NC_123456.1|', 3='gi|123456|ref|NC_123456.1| Strain name'
         -t <STR>   : the name of target genome, have to be consistent with (-n); default [NC_000962]
         -f <STR>   : (optional), filename, write the position in in db

Example:
ext_homoseq.pl -o subseq -d db.fa -n id input.fa
\n");
}

### 
__END__

Change log:

1. 2015-08-01
    Create this script

2. 2015-09-10
    Add a parameter: -t, move the target genome to the first place



### END OF FILE ###
