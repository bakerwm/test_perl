#!/usr/bin/perl -w

## parase the output of fastQC




##

use strict;
use warnings;
use File::Basename qw(basename);

## input: dir, contain zip files
my $inDir = shift or die "Need input a directory containning ZIP files:\n";
my @zips  = glob"$inDir/*.zip";

for my $z (@zips) {
    print basename($z) . "\t";
    my @info = parse_fastqc($z);
    print join("\t", @info) . "\n";
}


sub parse_fastqc {
    my $in = $_[0]; # *.zip file
    ## create a temp_dir
    my $tmp=`echo $in | md5sum | cut -d\" \" -f 1`;
    chomp($tmp);
    if( ! -d $tmp ){
        system "mkdir -p $tmp";
    }
    system "unzip -q -o -d $tmp $in";
    ##
    my $in_name = basename($in); $in_name =~ s/\.zip$//;
#    my $data = glob"$tmp/$in_name/fastqc_data.txt";
    my $data_file = `find $tmp -name "fastqc_data.txt"`;
    chomp($data_file);
    my @out = ('-', '-', '-', '-', '-', '-');
    if(-f $data_file) {
        @out = (); # clear objects in @out
        my %db = %{read_data($data_file)};
        push @out, get_basicstat($db{'Basic Statistics'});
        push @out, get_overrepresented($db{'Overrepresented sequences'});
    }
    ##
    system "rm -fr $tmp";
    ##
    return(@out);
}

sub read_data {
    my $in = $_[0]; # fastqc_data.txt
    $/ = ">>";
    open my $f_in, "< $in" or die "Cannot open file $in, $!\n";
    my @blocks = <$f_in>;
    close $f_in;
    $/ = "\n";
    ## parse each block
    my %db = ();
    for my $b (@blocks) {
        $b =~ s/\>\>$//;
        my @lines = split "\n", $b;
        my $head_line = shift(@lines);
        my $head_title = (split /\t+/, $head_line)[0];
        my $body = join "\n", @lines;
        $db{$head_title} = $body;
    }
    return(\%db);
}

sub get_basicstat {
    my $in = $_[0]; # string, Basic Statistics
    my @lines = split "\n", $in;
    my %out = ();
    for my $t (@lines) {
        my ($name, $content) = split "\t", $t;
        $out{$name} = $content;
    }
    ## Filename
    my @select = ();
    push @select, $out{'Filename'};
    push @select, $out{'Total Sequences'};
    push @select, $out{'Sequence length'};
    push @select, $out{'%GC'};
    return(@select);
}

sub get_overrepresented {
    my $in = $_[0]; # string, Overrepresented sequences
    my @out = ('-', '-');
    if((split/\t/,$in) > 1) {
        @out = ();
        my @tank_ad = ();
        my @tank_wt = ();
        my @lines = split "\n", $in;
        for my $t (@lines) {
            next if($t =~ /^\#/);
            my @tabs = split /\t/, $t;
            if($tabs[3] eq 'No Hit') {
                push @tank_wt, $tabs[0];
                push @tank_wt, $tabs[3];
            } else {
                push @tank_ad, $tabs[0];
                push @tank_ad, $tabs[3];
            }
        }
        ## which one to output
        @out = @tank_wt[0, 1];
        if(@tank_ad > 1) {
            @out = @tank_ad[0, 1];
        }
    }
    return(@out);
}

## EOF
