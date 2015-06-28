### !---------------------------
### replaced by chk_feature2count.pl

#!/bin/bash/perl -w
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

my $usage = "perl $0 -i infile -c 7 -n 1000000 > outfile";
my ($infile, $num, $help);
my $col = '0';
GetOptions(
    'input|i=s' => \$infile,
    'col|c=i' => \$col,
    'num|n=i' => \$num,
    'help|h' => \$help
    );

pod2usage(1) if($help);
die("$usage \n[Need input -i and -n]") unless (defined $infile && defined $num);
#die("$usage\n") if(@ARGV == 0);

## Parsing the input file
open F, "< $infile" or die "$!";

while(<F>){
    chomp;
    my @tabs = split "\t";
    die "[$.: input: $infile] is not a tab-separated file" if(scalar(@tabs) < 2);
    my $Ncol = $col - 1;
    my $count = $tabs[$Ncol];
    die "[$.: col: $col = $count] is not an integer" unless(is_integer($count));
    my $TPM = sprintf "%.2f", ($count * 1e6) / $num;
    print join "\t", (@tabs, $TPM), "\n";
}

close F;
#close OUT;

sub is_integer{
    my $t = shift(@_);
    my $flag = ($t=~/^[\+\-]?\d+$/)?"1":"0";
    return $flag;
}

__END__

=head1 NAME

c<count_to_TPM.pl> - Calculate the TPM value.

=head1 SYNOPSIS

# perl count_to_TPM.pl -i input.txt -c 7 -n 1000000 > out.txt

=head1 OPTIONS

=over 8

=item B<-i> 

A Tab-separated file with count number.

=item B<-c>

The <N> is the column that count number is in. last column [default: 0]

=item B<-n>

The total number of reads in this library.

=item B<--help>

Print this help

=head1 AUTHOR

Wang Ming <wangmcas@gmail.com>
