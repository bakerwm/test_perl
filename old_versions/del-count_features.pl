### !----------------------------
### bedtools counting is not right.(redundency)

#!/usr/bin/perl -w
use strict;
use warnings;

use Getopt::Std;
use File::Basename qw(basename dirname);
use File::Spec::Functions qw(catfile catdir);
use POSIX qw(strftime);

my %opts = ();
getopt("g:b:o:m:", \%opts);
die "Usage: perl $0 [-b] <Alignment.dir> [-o] <result.dir> [-m] <c|s>  <feature.dir>\n" if(@ARGV != 1);

## Tools
my $bedtools = '/home/wangming/software/bedtools2-2.21.0/bin/bedtools';
my $sort2bed = '/home/wangming/work/bin/sort2bed.pl';

## parse feature files
my $mRNA_bed;
my $AS_bed;
my $IGR_bed;
my $rnt_bed;

## parse BAM files
my @bam_files = <$opts{b}\/*.s.bam>;
die "check bam_files in $opts{b}" if(@bam_files < 1);
my $bam_list = join " ", (sort @bam_files);

mkdir $opts{o} unless -d $opts{o};
## parse feature files: mRNA, AS, rRNA, tRNA, IGR,
my $f_dir = shift;
my @f_files = <$f_dir\/*.bed>;

my $result = catfile($opts{o}, "genome_features.stat");
open my $res, "> $result" or die "$!";

print '[', &show_date, '] Start stat genome features', "\n";
foreach my $f (@f_files) {
    my $stat = '';
    if($f =~ /(mRNA|AS|rRNA|tRNA|IGR)\.bed$/){
        $stat = &run_bedtools($f);
    }else{
        next;
    }
    print $res $stat, "\n";
}
close $res;
print '[', &show_date, '] Finish.', "\n";

## 
sub run_bedtools {
    my $in = shift(@_);
    my ($type) = $in =~ /(\w+)\.bed/;
    my $count_out = catfile($opts{o}, "$type\_count.txt");
    my $cmd = join " ", ('bedtools', 'multicov', '-bams', $bam_list, '-bed', $in, '-s >', $count_out);

    if($opts{m} eq 'c') {
        print '[bedtools command:] ', $cmd, "\n";
        print '[', &show_date, '] Run bedtools - ', basename($in), "\n";
        system "$cmd &";
        return $cmd;    
    } elsif($opts{m} eq 's') {
        # stat sum
        print '[', &show_date, '] sum counts', "\n";
        my @sum_counts = &sum_count($count_out);
        # output
        my $line = join " ", ($type, @sum_counts);
        return $line;
    } else{
    
    }
}
##
sub sum_count {
    my $in = shift(@_);
    my @sum = ();
    open my $fh, "< $in" or die "$!";
    while(<$fh>) {
        chomp;
        my @tabs = split /\s+/;
        for(my $i = 0; $i < @bam_files; $i++){
            $sum[$i] += $tabs[$i+6];
        }
#        $sum[0] += $tabs[9];
#        $sum[1] += $tabs[10];
#        $sum[2] += $tabs[11];
#        $sum[3] += $tabs[12];
    }
    close $fh;
    return @sum;
}
##
sub show_date {
    my $date = strftime "%Y-%m-%d %H:%M:%S", localtime; 
    return $date;
}

