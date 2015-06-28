### !-------------------------------------------
### replaced by: samtools flagstat 

#!/usr/bin/perl -w
use strict;
use warnings;

use lib qw(/home/wangming/localperl/share/perl5);

use File::Basename qw(dirname basename);
use File::Spec::Functions qw(catdir catfile);
use File::Which qw(which);
use Cwd qw(cwd abs_path);
use Getopt::Long;
use Pod::Usage;
use POSIX qw(strftime);
use Data::Dumper;


my $bam_dir = '';
my $outstat = '';
my $outlog  = '';
my $samtools = '';
my $help = '';

GetOptions('bamdir|b=s' => \$bam_dir,
           'output|o=s' => \$outstat,
           'logout|l=s' => \$outlog,
           'samtools|s=s' => \$samtools,
           'help|h' => \$help
        );
pod2usage(1) if $help;

# Example:
# samtools flagstat in.s.bam > in.s.bam.stat
# Check samtools
my $tools = which('samtools');
if($samtools eq ''){
    if(defined $tools){
        $samtools = $tools;
    }else{
        pod2usage(-message=>"$0: Need input the path of samtools. [-s|--samtools path_to_samtools]") unless(defined $tools);
    }
}
die "File not fount: [$samtools]" unless -e $samtools;

# Check bam files
my @in_bams = <$bam_dir/*.s.bam>;
pod2usage(-message=>"$0: No *.s.bam files found under [-b $bam_dir].") if (@in_bams < 1);

# Start pasrsing BAM files
my %result = ();
my @log_info = ();
foreach my $b (@in_bams){
     my $b_name = basename($b);
     die "Not a sorted BAM file: [$b]" unless($b_name =~ /\.s\.bam/);
     my ($file_name) = $b_name =~ /(.*)\.s\.bam/;
     my $log = strftime "%Y-%m-%d %H:%M:%S, Statistic bam file: $b ...", localtime;
     push @log_info, $log;
#     print $log_info, "\n";
     my @tmp = `$samtools flagstat $b`;
#     die "The output of $b is not complete." if (@tmp != 13);
     my @numbers = &parse_num(@tmp);
     my $percent = sprintf"%.4f",($numbers[8]/$numbers[0]);
     $result{$file_name} = join" \t", (@numbers[0,8], $percent);
}

$outstat = 'out.stat' if ($outstat eq '');
$outlog = 'out.log' if ($outlog eq '');

open F1, "> $outlog" or die "$outlog, $!";
open F2, "> $outstat" or die "$outstat, $!";

print F1 join"\n", @log_info, "\n";
#my $stat_header = join"\t", ('Name', 'Total', 'Total-NotQC', '2nd-1', '2nd-2',
#                             'sup-1', 'sup-2', 'dup-1', 'dup-2', 'Map-1', 'Map-2',
#                             'pair-1', 'pair-2', 'read1-1', 'read1-2',
#                             'read2-1', 'read2-2', 'prop-1', 'prop-2', 'Mp-1', 'Mp-2',
#                             'sing-1', 'sing-2', 'dchr-1', 'dchr-2', 'dchr-3', 'dchr-4'
#                             );
my $stat_header = join"\t", ('Name', 'Total', 'Mapped', 'Ratio');
print F2 $stat_header,"\n";                             
foreach my $f (sort keys %result){
    print F2 join "\t", ($f, $result{$f}), "\n";
}

close F1;
close F2;

#print Dumper(%result);

sub parse_num{
    my @ta = ();
    foreach(@_){
        my ($num1, $num2) = /(^\d+)\s\+\s(\d+)\s/;
        push @ta, ($num1, $num2);
    }
    return @ta;
}



__END__

=head1 NAME

c<statMapping.pl> - Statistic the mapping of input reads.

=head1 SYNOPSIS

# perl statMapping.pl -b bamfile_dir -o out.stat -l out.log

=head1 OPTIONS

=over 8

=item B<-b>, str, B<--bamdir>

The directory that contain BAM files, with the following naming
criteria:[Name].s.bam

=item B<-o>, str, B<--output>

The file to store the statment of BAM files. [default=out.stat]

=item B<-l>, str, B<--logfile>

Save the details of what this script done to a log file.
[default=out.log]

=item B<-s>, str, B<--samtools>

The path of samtools. This script will search the PATH for samtools first, and it will
reminde you to input the path of samtools if it does not found it.

=item B<-h>, B<--help>

Show this help.

=back

=head1 AUTHOR

Wang Ming <wangmcas@gmail.com>
Nov-29, 2014

=cut
