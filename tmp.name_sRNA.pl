#!/usr/bin/env perl

##################################################
# Rename the ncRNAs according to the publication:
#
# Lamichhane G, et al., Definition and annotation 
# of (myco)bacterial non-coding RNA, Tuberculosis 
# (2012), 
# http://dx.doi.org/10.1016/j.tube.2012.11.010
##################################################

use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;

name_sRNA();
#exit(1);

sub name_sRNA {
    my %opts = (n => 'sRNA');
    getopts("a:l:n:o:", \%opts);
    usage() if(@ARGV == 0 && -t STDIN);
    die("[-a] Need *.ptt file (NCBI)\n") if(! defined $opts{a});
    # read input txt
    my @lines = ();
    while(<>) {
        chomp;
        next if(/^\#|^\s*$/);
        push @lines, $_;
        die("Input should contain at least 12 columns\n") if((split /\t/,$_) < 12);
    }
    # type
    my $type = '';
    if($opts{n} =~ /^sRNA$/) {
        $type = 'nc';
    }elsif($opts{n} =~ /^UTR$/) {
        $type = 'utr';
    }elsif($opts{n} =~ /^IM$/) {
        $type = 'im';
    }elsif($opts{n} =~ /^OTR$/) {
        $type = 'otr';
    }elsif($opts{n} =~ /^del$/) {
        $type = 'del'
    }else {
        die("[-n $opts{n}] Unknown type, expect: sRNA, UTR, IM\n");
    }
    # rename
    my @new_file = locate_tags(\@lines, $type, $opts{a});
    # ouput
    my $fh_out = *STDOUT;
    if(defined $opts{o}) {
        open $fh_out, "> $opts{o}" or die "Cannot write to $opts{o}, $!\n";
    }
    print $fh_out join("\n", @new_file) . "\n";
}
 
sub locate_tags {
    my @list = @{$_[0]};
    my $type = $_[1]; # nc, utr, im
    my $ptt  = $_[2];
    # split into two categories: IGR + AS
    my %hit = ();
    for my $i (@list) {
        my @tabs = split /\t/, $i, 13;
        if($tabs[11] eq 'IGR') {
            push @{$hit{'IGR'}->{$tabs[6]}}, $i;
        }elsif($tabs[11] =~ /AS|PM|IM/) {
            my $g = $tabs[6];
#            if($tabs[6] ne $tabs[8] && $tabs[7] > 0 && $tabs[9] < 0) {
            if($tabs[6] ne $tabs[8] && $tabs[9] < $tabs[7]) {
                $g = $tabs[8];
                if($tabs[7] < 0 && $tabs[9] < 0) {
                    $g = $tabs[6];
                }
            }
            push @{$hit{'AS'}->{$g}}, $i;
        }else {
            push @{$hit{'NA'}->{$tabs[6]}}, $i;
#push @{$hit{'NA'}}, $i;
        }
    }
    # assign new id to each seq
    # ncRv, ncRv1, A->..., c, utr
    # IGR
    my @reports = ();
    for my $k (sort keys %{$hit{'IGR'}}) {
        next if($k eq '');
        my @subs = sort_pos( @{$hit{'IGR'}->{$k}} );
        my $count = 1;
        for my $s (@subs) {
           my @das = split /\t/, $s;
           my $old_id = $das[0];
           $das[0] = new_name($k, $type, 'IGR', $count, $das[5], $ptt);
           push @reports, join("\t", @das, $old_id);
           $count ++;
        }
    }
    # AS
    for my $j (sort keys %{$hit{'AS'}}) {
        my @subas = sort_pos( @{$hit{'AS'}->{$j}} );
        my $count = 1;
        for my $p (@subas) {
            my @das = split /\t/, $p;
            my $old_id = $das[0];
            $das[0] = new_name($j, $type, 'AS', $count, $das[5], $ptt);
            push @reports, join("\t", @das, $old_id);
            $count ++;
        }
    }
    # NA
    for my $n (sort keys %{$hit{'NA'}}) {
        my @subna = sort_pos( @{$hit{'NA'}->{$n}} );
        my $count = 1;
        for my $m (@subna) {
            my @das = split /\t/, $m;
            my $old_id = $das[0];
            $das[0] = new_name($n, $type, 'OTR', $count, $das[5], $ptt);
            push @reports, join("\t", @das, $old_id);
            $count ++;
        }
    }
    return(@reports);
}

sub sort_pos {
# sort input seqs by the position (col-4)
    my @out = ();
    my %t   = ();
    for my $i (@_) {
        my @tabs = split /\t/, $i;
        $t{$tabs[3]}->{$tabs[4]}->{$tabs[5]} = $i;
    }
    for my $j (sort {$a<=>$b} keys %t) {
        for my $k (sort {$a<=>$b} keys %{$t{$j}}) {
            for my $l (sort keys %{$t{$j}->{$k}}) { # +/-
                push @out, $t{$j}->{$k}->{$l};
            }
        }
    }
    return(@out);
}

sub guess_prefix { # guess the prefix of gene name: Rv0001 "Rv", MRA_0001 "MRA", Rvnrn1 "Rv"
    my $in = $_[0]; # ptt file
    my %h  = ();
    open my $fh_in, "< $in" or die "Cannot opne $in, $!\n";
    my $count = 1;
    while(<$fh_in>) {
        chomp;
        last if($count > 200);
        next if(/^\#|^\s*$/);
        my @tabs = split /\t/, $_;
        next if(@tabs < 9);
        # estimate the length of prefix
        if($tabs[5] =~ /\_/) {
            my $p1 = (split /\_/, $tabs[5])[0];
            my $p1_length = length($p1);
            next if($p1_length eq '');
            $h{$p1_length} ++;
        }else {
            my ($p2) = $tabs[5] =~ /([a-zA-Z]+)\d+/;
            my $p2_length = length($p2);
            next if(! defined $p2_length);
            $h{$p2_length} ++;
        }
        $count ++;
    }
    close $fh_in;
    my $main_length = 1;
    my $max_count   = 1;
    for my $i (sort {$a<=>$b} keys %h) {
        if($h{$i} > $max_count) {
            $main_length = $i;
            $max_count   = $h{$i};
        }
    }
    return $main_length;
}

sub new_name {
    my $id     = $_[0];
    my $type   = $_[1]; # nc, utr, im
    my $cat    = $_[2]; # IGR or AS
    my $order  = $_[3];
    my $strand = $_[4];
    my $ptt    = $_[5];
    #
    my @suffix = ('A'..'Z', 'AA'..'ZZ'); # 676 elements
    my $id_prefix_width = guess_prefix($ptt);
    $id =~ s/c$//;
    my ($p1, $p2);
    if($id =~ /\_/) {
        ($p1, $p2) = split /\_/, $id;
    }else {
        ($p1, $p2) = $id =~ /([a-zA-Z]{$id_prefix_width})(\w+)/;
    }
    # each parts
    my $n1 = $type; # nc, utr, im
    my $n2 = $p1; # Rv, MRA
    my $n3 = ''; # 1 or ''
    my $n4 = $p2; # 0001
    my $n5 = ''; # rank, A-Z
    my $n6 = ''; # c or ''
    # 
    $n3 = '1' if($cat eq 'IGR');
    $n5 = $suffix[$order - 1]; # A-Z...
    $n6 = 'c' if($strand eq '-');
    # ncRv10001Ac
    return join('', $n1, $n2, $n3, $n4, $n5, $n6);
}

sub usage {
    die("
Usage: name_sRNA.pl [options] <in.txt|STDIN>

Options: -a <STR>           : PTT annotation file (NCBI)
         -n <STR>           : name of the sRNA (sRNA|UTR|IM|OTR|del) [sRNA]
         -o <STR>           : (optional) redirect results to file
         
Example:
nam_sRNA.pl -a NC_000962.ptt -n sRNA in.txt > out.txt            
\n");
}

sub parse_ptt { # save all the genes in the order of corrdinations.
    my $in = $_[0];
    my %info = ();
    my %check = ();
    open my $fh_in, "< $in" or die "Cannot open $in, $!\n";
    while(<$fh_in>) {
        chomp;
        next if(/^\#|^\s*$/);
        next unless(/^\d+\.\.\d+\t/);
        my ($pos, $strand, $id) = (split /\t/)[0,1,5];
        my ($start, $end) = split /\.+/, $pos;
        if(exists $check{$pos}) {
            die("multiple features in the same genomic position:\n$check{$pos}\n$_\n");
        }
        $check{$pos}  = $_;
        $info{$start} = $id;
    }
    close $fh_in;
    return(\%info);
}


__END__

Change log:

2015-07-06
  v0.1    begin script
          

  Need to do:
  PTT: only used for evalulate the name of cds. (further plan: locate sRNA)

  v0.2  
  warning: this script will discard the duplicate seqs (sort by: start,end,strand)
