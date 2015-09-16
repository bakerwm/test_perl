#!/usr/bin/env perl

#################################################
# tmp.UTRtag.pl 
#
# add utr tags to txt file,
# add genes to txt file (sRNA associated)
#
#################################################

use strict;
use warnings;
use Getopt::Std;

utr_tag();
exit(1);

###
sub utr_tag {
    my %opts = ();
    getopts("o:", \%opts);
    
    usage() if(@ARGV == 0 && -t STDIN);
    my $fh_out = *STDOUT;
    if(defined $opts{o}) {
        open $fh_out, "> $opts{o}" or die "Cannot write to $opts{o}, $!\n";
    }
    ###
    while(<>) {
        chomp;
        my ($gene, $tag) = name_utr($_);
        print $fh_out $_ . "\t" . $gene . "\t" . $tag . "\n";    
    }

}
###
# AS, IGR, PM
# UTR/ 5', 3'
sub name_utr {
    my $in = $_[0];
    my @tabs = split /\t/, $in;
    die("[$in] less than 12 columns\n") if( @tabs < 12 );
    my ($g1, $d1, $g2, $d2, $dir, $des) = @tabs[6..11];
    my ($s1, $s2, $s3) = (split /\//, $dir)[1..3];

    my $tag  = ''; # or utr5
    my $gene = '';
    if($s1 eq $s2) {
        if($s2 eq $s3) {
            if($d1 < $d2) {
                $gene = $g1;
                $tag  = ($s2 eq '+')?'utr3':'utr5';
            }else {
                $gene = $g2;
                $tag  = ($s2 eq '+')?'utr5':'utr3';
            }
            if($des =~ /^PM\d+$/) {
                $tag  = 'utr35';
                $gene = $g1 . ',' . $g2; 
            }
        }else {
            $gene = $g1;
            $tag  = ($s2 eq '+')?'utr3':'utr5';
        }
    }else {
        if($s2 eq $s3) {
            $gene = $g2;
            $tag = ($s2 eq '+')?'utr5':'utr3';        
        }else {
            $tag  = '-';
            $gene = '-';
        }
    }
    return ($gene, $tag);
}

sub usage {
    die("
Usage: UTRtag.pl [options] <in.txt|STDIN>;
            
Options: -o     : redirect results to file (optional)
         
<in.txt>
have to contain the 12-col part: (strand and gene position)

Add [gene]+ [tag] to the tail of each line
\n")
}


