#!/usr/bin/perl  -w
use strict;
use warnings;
use Data::Dumper;
####################################
# Date: Feb 28, 2013
####################################

my $mapping_dir = shift || die "Need input mapping_dir: $! \n";
my $input_file  = shift || die "Need input sRNA text file: $! \n";
my @mapping_files = glob("$mapping_dir*mapping.txt");
#die "Usage: perl  Cal_sRNA_end.pl  mapping_dir/  input.txt \n" if(@ARGV < 2);
die "Check $mapping_dir for *_mapping.txt files\n" if(@mapping_files < 1);

# Read and store mapping files.
my %h;
foreach my $f (@mapping_files){
    $f =~ /\/(\d+)\w+\_(\w+)\.mapping\.txt/; # /Coverage/140PE_3end.mapping.txt
    my ($lib, $type) = ($1, $2);
  # Read mapping file
    open F,$f || die "Cannot open file $f: $! \n";
    while(<F>){
        chomp;
        my ($position, $coverage, $str) = split(/\t/,$_);
        my $key = join"\:",($type, $str);
        $h{$key}->{$lib}->{$position} = $coverage;
    }
    close F;
}

# Check Cov_* directories.
unlink "Cov_*/*.txt" if (-e "Cov_*/*.txt");
rmdir  "Cov_*/" if (-e "Cov_*");

open F, $input_file || die "Cannot open file $input_file: $! \n";
while(<F>){
    chomp;
    my @ps = split(/\t/,$_);
    &Report($.) if(@ps < 6);
    my $start = int($ps[3]/100) * 100;
    my $end   = int($ps[4]/100 + 1) * 100;
    my ($ID, $SeqStrand) = ($ps[0], $ps[5]);
    $SeqStrand =~ s/\r//;

    foreach my $i(sort keys %h){
        $i =~ /(\w+)\:(.)/;
        my ($type, $In_strand) = ($1, $2);
        next if($In_strand ne $SeqStrand);
        mkdir ("Cov_$type", 0755);
        my $outfile = "Cov_$type\/$ID\_$type\.txt";        
        open OUT,">$outfile" || warn "Cannot open file $outfile: $! \n";
        for(my $k=$start; $k<=$end; $k++){
            my @hit_cov = ();
            foreach my $m(sort {$a<=>$b} keys %{$h{$i}}){
                my $hit = (exists $h{$i}->{$m}->{$k})?$h{$i}->{$m}->{$k}:'0';
                push @hit_cov,$hit;
            }
            my $hit_out = join"\t",($k,@hit_cov);
            print OUT $hit_out,"\n";
        }
        close OUT;    
    }
}
close F;

sub Report{
    my $line = shift(@_);
    print "Check line $line in $input_file, which should be 6-line and +/- at line-5 \n";
    exit(1);
}
