### !------------------------------
### replaced by: chk_parseBAM.pl

#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename qw(basename dirname);
use File::Spec::Functions qw(catfile catdir);


my $list = shift or die "Usage: combine_TPMs.pl <TPM.list>  \n";

open IN, "<$list" or die "$!";
while(<IN>){
    chomp;
    my $in = $_;
    die "[$in] not exist\n" unless -e $in;
    my $tpm = $in;
    my $txt = $in;
    $txt =~ s/\_TPM\.bed/.txt/;
    my $rnaz = catfile(catdir(dirname($_), "RNAz_out"), "best_RNAz.bed");
    die "[$rnaz] not found. Need execute: sh run_rnaz.sh\n" unless -e $rnaz;

    my $out_name = basename($in);
    $out_name =~ s/\_TPM\.bed/_combine.txt/;
    my $out = catfile(dirname($in), $out_name);
    
    print join"\n", ('Input files:', $in, $txt, $rnaz, $out),'';
    combine_files_by_id($txt, $tpm, $rnaz, $out);

}
close IN;

sub combine_files_by_id {
    my ($txt, $tpm, $rnaz, $out) = (@_);
    die "Need [4] input files" unless(@_ == 4);
    my %sort = ();
    my %tpm = ();
    my %rnaz = ();
    %sort = read_txt($txt);
    %tpm  = read_bed($tpm);
    %rnaz = read_rnaz($rnaz);

    open OUT, ">$out" or die "$!";
    my $count_tpm = '-'. "\t". '-';
    foreach my $id (keys %sort) {
        if(exists $tpm{$id}) {
            $count_tpm = $tpm{$id};
        }else{
            die "[$id] not found count_tpm: [$tpm, $txt]\n";
        }
        my $z_score = '-';
        $z_score = $rnaz{$id} if(exists $rnaz{$id});
        print OUT join "\t", ($id, $sort{$id}, $count_tpm, $z_score);
        print OUT "\n";
    }
print "Write output:\n $out\n\n";
}

sub read_txt {
    my $in = shift(@_);
    my %sort = ();
    open F, "<$in" or die "$!";
    while(<F>) {
        chomp;
        next if(/^\s*$/);
        my @tabs = split /\t/, $_;
        my $id = shift(@tabs);
        my $id_new = $id;
        $id_new = (split /\,/, $id)[0] if($id =~ /\,/);
        $sort{$id_new} = join "\t", (@tabs, $id);
    }
    close F;
    return %sort;
}

sub read_bed {
    my $in = shift(@_);
    my %bed = ();
    open F, "<$in" or die "$!";
    while(<F>) {
        chomp;
        next if(/^\s*$/);
        my ($id, $count, $tpm) = (split /\t/, $_)[3,6,7];
        $id = (split /\,/, $id)[0] if($id =~ /\,/);
        $bed{$id} = $count. "\t". $tpm;
    }
    close F;
    return %bed;
}

sub read_rnaz {
    my $in = shift(@_);
    my %rnaz = ();
    open F, "<$in" or die "$!";
    while(<F>) {
        chomp;
        next if(/^\s*$/); # skip blank lines
        my ($id, $score) = (split /\t/, $_)[0,4];
        my $id_new;
        if($id =~ /LIB/){
            $id_new = (split /\_LIB/, $id)[1];
            $id_new = 'LIB'.$id_new;
        }elsif($id =~ /Seed/) {
            $id_new = (split /\Seed/, $id)[1];
            $id_new = 'Seed'.$id_new;
        }else{
        }
        $rnaz{$id_new} = $score;
    }
    close F;
    return %rnaz;
}

