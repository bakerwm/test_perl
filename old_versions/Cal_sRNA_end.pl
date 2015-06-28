### !---------------------------------------------------
### delete


#!/usr/bin/perl -w
use strict;
use Data::Dumper;

sub help{
    print STDERR <<EOF;

Usage: perl  Cal_End_cov.pl  mapping_dir  infile
Note:
1. Mappding_dir include files: 45SE_3end.mapping.txt ...
2. Input files should be in SORT format:
ID*  Exp Length  Begin*   End* Strand*   (* is required)

EOF
}

my $mapping_dir = shift or die &help;
my $infile      = shift or die &help;
my @map_files   = glob"$mapping_dir*mapping.txt";

if(@map_files == 0){
    print "Warn: $mapping_dir does not contain *mapping files.\n";
    &help();
    exit(1);
}

my %mapping;
foreach my $map (@map_files){
    $map =~ m/(\d+[P,S]E)\_(\w+)\.mapping\.txt/;    # eg: 45SE_3end.mapping.txt
    my ($lib, $type) = ($1, $2);

    open F, "<$map" or die;
    while(<F>){
        chomp;
        my($position, $cov, $strand) = split (/\t/);
        $mapping{$type}->{$lib}->{$strand}->{$position} = $cov;
    }
    close F;
}

# Remove directories.
unlink "Cov_*/*";
rmdir  "Cov_*/";

my $readme;
# Calculate coverage.
open F,  "< $infile"    or die;
open Chk,"> Readme.txt" or die;
while(<F>){
    chomp;
    my @ps     = split(/\t/);
    my $id     = $ps[0];
    my $start  = int($ps[3]/100) * 100;
    my $end    = (int($ps[4]/100) + 1) * 100;
    my $strand = $ps[5];
    foreach my $In_type (sort keys %mapping){   # 3end or 5 end or Poscov
        mkdir ("Cov_$In_type", 0755);
        open OUT,"> Cov_$In_type/$id\_$In_type\.txt" or die;
        
        my %pos_cov;
        foreach my $In_lib(sort keys %{$mapping{$In_type}}){    # 45SE, 81SE, ...
            $In_lib =~ m/(\d+)[S,P]E/;
            my $tag =  $1; # 45, 81, 140, 200
            foreach my $In_strand (sort keys %{$mapping{$In_type}->{$In_lib}}){ # strand ...
                next unless($In_strand eq $strand);
                for(my $i = $start; $i <= $end; $i++){
#                    if(exists $mapping{$In_type}->{$In_lib}->{$In_strand}->{$i}){
#                        $pos_cov{$tag}->{$i} = $mapping{$In_type}->{$In_lib}->{$In_strand}->{$i};
#                    }else{
#                        $pos_cov{$tag}->{$i} = 0;
#                    }
                    $pos_cov{$tag}->{$i} = (exists $mapping{$In_type}->{$In_lib}->{$In_strand}->{$i})?$mapping{$In_type}->{$In_lib}->{$In_strand}->{$i}:0;
                }
            
            }
        }

        for(my $i = $start; $i <= $end; $i++){
            my @check = ();
            my @out   = ();
            push @out,($i);
            push @check,($i);
            foreach my $g (sort{$a<=>$b} keys %pos_cov){
                push @out,($pos_cov{$g}->{$i});
                push @check,($g);
            }
            my $line = join"\t", @out;
               $readme  = join"\t", @check;

            print OUT $line,"\n";
#            print $chk,"\n";
        }
        close OUT;
    }
}
close F;

$readme = "The headline for each files in Cov_*/ is:
Position    lib_1   lib_2   lib_3   lib-4\n".$readme;

print Chk $readme,"\n";
