#!/usr/bin/perl

use strict;
use warnings;

use LWP::UserAgent;
#my $ua = LWP::UserAgent->new;
#my $req = HTTP::Request->new(GET => 'http://www.ebi.ac.uk/QuickGO/GAnnotation?protein=P12345,Q4VCS5&format=tsv&col=proteinID,proteinSymbol,evidence,goID,goName,aspect,ref,with,from');
#my $res = $ua->request($req, "annotation.tsv");

#open (FILE, 'annotation.tsv');
#my $head = <FILE>;
#while (<FILE>) {
#    chomp;
    my ($proteinID, $proteinSymbol, $evidence, $goID, $goName, $aspect, $ref, $with, $from) = split//, '-' x 9;
#    ($proteinID, $proteinSymbol, $evidence, $goID, $goName, $aspect, $ref, $with, $from) = split(/\t/);
#    print "$proteinID => $goID ($aspect: $goName) $evidence  $ref  $with  $from\n";
#}
#close FILE;
##
#exit;

    print join("\t", $proteinID, $proteinSymbol, $from) . "\n";
