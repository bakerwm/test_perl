### !----------------------------------------------
### replaced by: chk_seq2rnaz.pl

#!/usr/bin/perl -w
use strict;

use lib qw(/home/wangming/localperl/share/perl5
           /home/wangming/localperl/lib/5.16.1
           /home/wangming/localperl/lib/site_perl/5.16.1
        );

use File::Basename qw(dirname basename);
use File::Spec::Functions qw(catdir catfile);
use Cwd qw(cwd abs_path);

sub usage{
    print << "EOF";
Usage: perl RNAzStat_parse.pl file1.stat file2.stat file3.stat ...

Parse the stat output of rnazBEDstats.pl.

EOF
exit(1);
}

&usage if(@ARGV < 1);
my @infiles = @ARGV;

my %result = ();
foreach my $f (@infiles){
    warn "File not found: $f" unless -e $f;
    my $f_name = basename($f);
    die "Input file name should be: [*.stat], $f " unless($f_name =~ s/\.stat$//g); 
    open F, $f or die;
    while(my $in = <F>){
        next if($in =~ /^\s+$/);
        chomp($in);
        $in =~ s/\s+//g;
        my ($item, $num) = split /\:/, $in;
        $result{$f_name}->{$item} = $num;  
    }
    close F;
}

my @out = ();
my $header = "";
foreach my $i (keys %result){
    my @header = keys(%{$result{$i}});
    my @num    = values(%{$result{$i}});
    push @out, join"\t", ("cat", @header);
    push @out, join"\t", ($i, @num);
}
print join"\n", @out;
print "\n";







