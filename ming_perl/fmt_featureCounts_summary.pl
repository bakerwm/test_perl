#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename qw(basename);

###########################################################################
# This script is designed to format the ouput of featureCounts, make the  #
# results easy to use.                                                    #
#                                                                         #
###########################################################################

sub usage {
    die("
Usage: $0 <summary_dir>
\n");
}

usage() if(@ARGV != 1);
my $in_dir  = $ARGV[0];
my @in_sums = glob("$in_dir\/*.summary");
die("[*.summary] files not found in: $in_dir\n") if(@in_sums == 0);

my %fmt = ();
my @features = ();
for my $f (sort @in_sums) {
    my ($f_name) = basename($f) =~ /(.*)(_count)?(\.txt)?\.summary/;
    $f_name =~ s/\_count$//;
    push @features, $f_name; ### feature names
    my @f_stats  = parse_stats($f);
    for my $fs (@f_stats) {
        my ($bam, $mapped, $unmapped) = split /\t/, $fs;
        $fmt{$bam}->{map}->{$f_name} = $mapped;
        $fmt{$bam}->{un}->{$f_name}  = $unmapped;
    }
}

my $header = join("\t", "bam_file", @features, @features);
print $header . "\n";

for my $b (sort keys %fmt) {
    my @nums = ($b);
    for my $f1 (sort keys %{$fmt{$b}->{map}} ) {
        push @nums, $fmt{$b}->{map}->{$f1};
    }
    for my $f2 (sort keys %{$fmt{$b}->{un}}) {
        push @nums, $fmt{$b}->{un}->{$f2};
    }
    print join("\t", @nums) . "\n";
}

sub parse_stats {
    my $in_file = $_[0];
    open my $fh_f, "< $in_file" or die "Cannot open file, $in_file, $!\n";
    my @lines = <$fh_f>;
    close $fh_f;
    my $line  = join("", @lines);
    ### split file by bam file 
    my @stats = ();
    my @segs  = split /Status/, $line;
    shift(@segs);
    for my $s (@segs) {
        my @frags = split /\n/, $s;
        my $bam_file = shift(@frags);
        $bam_file    =~ s/\s+//g;
#        $bam_file    = basename($bam_file); ### optional
        my $mapped   = (split /\s+/, shift(@frags))[1];
        my $unmapped = 0;
        for my $fg (@frags) {
            my ($type, $num) = (split /\s+/, $fg)[0, 1];
            next if($type eq 'Unassigned_Unmapped');
            $unmapped += $num;
        }
        push @stats, join("\t", $bam_file, $mapped, $unmapped);
    }
    return(@stats);
}

### END OF FILE ###
