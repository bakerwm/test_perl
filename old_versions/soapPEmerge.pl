#!/usr/bin/perl  -w
use strict;

sub help{
        print STDERR <<EOF;
Usage: perl soapPE2txt.pl  <*PE.soap>  <out.txt> 
Description:
This script is designed to merge PE soap two reads into one.
Options:
    *PE.soap        : soap PE result.
    out.txt         : merge result.
    delete.txt      : delete merge results that too long (>300) or 
                      too short (<50). Default: "*PE.soap.delete"
    repeat.txt      : The reads that have repeat hits.
                      "*PE.soap.repeat"
EOF
}
my ($PE, $out, $delete, $repeat);
$PE     = shift or die &help; # soap PE result.
$out    = shift or die &help;
$delete = shift or $delete = "$PE\.delete";
$repeat = shift or $repeat = "$PE\.repeat";

open F,  "< $PE"     or die;
my %h;
while(<F>){
    chomp;
    my @ps   = split(/\t/,$_);
    my $line = join"\,",@ps;
#    my $line = join"\,",($ps[5], $ps[6], $ps[8]); # length, strand, genome_id, begin
    push @{$h{$ps[0]}}, ($line);
}
close F;

open OUT,"> $out"    or die;
open DEL,"> $delete" or die;
open REP,"> $repeat" or die;
#print DEL "read_2\tPair_len\tread_1_len\tread_2_len\n";
foreach(keys %h){
    my @check = @{$h{$_}};
    next if($_ =~ /\/1$/); # skip read_1
    if(@check > 1){ # skip the reads have repeat hits.
        print REP "$_\n"x@check;
        next;
    }
    
    my $r_1 = my $r_2 = $_;
    my @ts  = split(/\,/, $h{$r_2}[0]);
       $r_1 =~ s/\/2$/\/1/g;
    my @ps  = split(/\,/, $h{$r_1}[0]);

    my ($read_1_length, $read_1_strand, $read_1_begin) = ($ps[5], $ps[6], $ps[8]);
    my ($read_2_length, $read_2_strand, $read_2_begin) = ($ts[5], $ts[6], $ts[8]);
    my  $read_1_end = $read_1_begin + $read_1_length - 1;
    my  $read_2_end = $read_2_begin + $read_2_length - 1;

    my ($pair_begin, $pair_end) = &min_max($read_1_begin, $read_1_end, $read_2_begin, $read_2_end);
    my  $pair_length = $pair_end - $pair_begin + 1;
    
    if($pair_length < 50 || $pair_length > 300){
        print DEL $r_2,"\t",$pair_length,"\t",$read_1_length,"\t",$read_2_length,"\n";
        next;
    }
    ($ts[5], $ts[8]) = ($pair_length, $pair_begin);
    my $p_out = join"\t", @ts;
#    my $p_out = join"\t",($r_2, "Exp", $pair_length, $pair_begin, $pair_end, $read_2_strand);
    print OUT $p_out,"\n";
}
close DEL;
close OUT;

sub min_max{
    my $min = shift;
    my $max = shift;
    foreach(@_){
        $min = ($min<$_)?$min:$_;
        $max = ($max>$_)?$max:$_;
    }
    return ($min, $max);
}
