#!/usr/bin/perl  -w
use warnings;
use strict;
use Getopt::Std;

sub help{
    print STDERR <<EOF
Usage: perl  get_sRNA.pl  -i  H37Rv_coverage.n  -c 0.5  -s +  -o  H37Rv_sRNA.n
EOF
}

# update: 2014-11-06, wangming

my %opts = ();
getopt("i:c:s:o:", \%opts);
if(!defined $opts{i} || !defined $opts{s} || !defined $opts{o}){
    &help();
    exit(1);
}
$opts{c} = (defined $opts{c})?$opts{c}:'0.5';
my $label = ($opts{s} eq '+')?'_p':'_n';
my ($min_cov, $min_length) = (100, 16);

open F, $opts{i} or die;
my %cov = ();
my %pos = ();
my $genome_length = 0;
foreach my $i(<F>){
    chomp($i);
    my @t = split(/\t/,$i);
    $cov{$t[2]}->{$t[1]} = 1;
    $pos{$t[1]} = $t[2];
    $genome_length = ($genome_length < $t[1])?$t[1]:$genome_length;
}
close F;

open OUT,"> $opts{o}" or die;
my %hit = ();
my $num = 1;
foreach my $j (sort {$b<=>$a} keys %cov){
    foreach my $k (sort {$a<=>$b} keys %{$cov{$j}}){
        my ($begin_pos, $begin_cov, $end_pos, $end_cov) = (1,1,1,1);
# Search toward 5' end.
        next if(exists $hit{$k});   
        my ($max_pos, $max_cov) = ($k, $j);
        next if($max_cov < $min_cov);
        for(my $m=$k;$m>=0;$m--){
            my ($check_pos, $check_cov) = ($m, $pos{$m});
            if(exists $hit{($check_pos - 1)}){
                ($begin_pos, $begin_cov) = ($m, $pos{$m});
                $hit{$m} = 1;
                last;
            }else{
                if($check_pos == 1 || $pos{($check_pos - 1)} < $opts{c} * $max_cov){   # check_cov < 1/2 of max_cov
                    ($begin_pos, $begin_cov) = ($m, $pos{$m});
                    $hit{$m} = 1;
                    last;
                }else{
                    $hit{$m} = 1;
                }
            }
        }
# Search toward 3' end.
        for(my $n=$k;$n<=$genome_length;$n++){
            my ($check_pos, $check_cov) = ($n, $pos{$n});
            if(exists $hit{($check_pos + 1)}){
                ($end_pos, $end_cov) = ($n, $pos{$n});
                $hit{$n} = 1;
                last;
            }else{
                if($check_pos == $genome_length || $pos{($check_pos + 1)} < $opts{c} * $max_cov){
                    ($end_pos, $end_cov) = ($n, $pos{$n});
                    $hit{$n} = 1;
                    last;
                }else{
                    $hit{$n} = 1;
                }
            }
        }
        my $number = sprintf"%04d",$num;
        my $seed_length = $end_pos - $begin_pos + 1;
        next if($seed_length < $min_length);
        my $out = join"\t",("Seed$number$label", $opts{s}, "$begin_pos\:$begin_cov", "$max_pos\:$max_cov", "$end_pos\:$end_cov", $seed_length);
        print OUT $out,"\n";
        $num ++;
    }
}
close OUT;
