#!/usr/bin/env perl

##############################################
# Convert the ptt to 
# id list.
#
# id gene Length Location Strand COG Product
# 
# Need *.ptt and *.gff in the same dir
##############################################

sub usage {
    die("
Usage: ptt2id.pl <in.gff|STDIN>

Extract gene info from GFF file 
\n");
}

usage() if(@ARGV == 0 && -t STDIN);

#while(<>) {
#    chomp;
#    next if(/^\#|^\s*$/);
#    next unless(/^\d+\.\.\d+\t/);
#    my @tabs = split /\t/;
#    my ($s, $e) = split /\.+/, $tabs[0];
#    my $length  = $e - $s + 1;
#    my $tag = ($tabs[4] eq '-')?$tabs[5]:$tabs[4];
#    print join("\t", $tag, $tabs[5], $length, $s, $e, $tabs[1], $tabs[3], $tabs[8]). "\n";
#}

my $feature = 'gene';

while(<>) {
     chomp;
     next if(/^\#|%\s*$/);
     next unless();

}






