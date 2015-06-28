### !-------------------------------
### replaced by: bedtools merge

#!/usr/bin/perl -w
use warnings;
use strict;
use Cwd qw/abs_path/;
use File::Basename qw/basename dirname/;
use Getopt::Std;
use Data::Dumper;

sub usage{
    print STDERR <<EOF;
Usage: perl merge4all.pl file1 file2 ...

Output: merge4all.txt

Note:
Merge RNA transcripts from different libraries that have at least one base overlap
according to the coordinations of each transcript.

Input file name:  *01_tag.txt
<Sample Name><Lib Name>_tag.txt 
<Lib Name> should be Integer[01-09], which represents: 
01=18-40nt, 02=40-48nt, 03=140nt, 04=200nt, ...

Input file format:
Tab-delimited, with at least 6 columns: (* are required)
<ID*> <Note> <Note> <Begin*> <End*> <Strand*>

Criteria for merging:
1. Have an overlap more than 10 nt or exceed the 40% of either of the transcript.
2. The reads count between less than 100-fold.

EOF
exit(1);
}

# Show the help
&usage if(@ARGV < 1);

################################################################################
# Create temp workdir: merge_tags
#my $out_dir = "04.Merge_tags";
#mkdir $out_dir unless -d $out_dir;
#my $out_dir_abs_path = abs_path($out_dir);

################################################################################
# Cat all files together
mkdir "tmp" unless -d "tmp";
my @tmps = <tmp\/tmp*\.txt>;
unlink @tmps; 

# Add "lib" name to the ID of each sequence
foreach my $f (@ARGV) {
    my ($f_lib) = $f =~ /(\d\d)\_tag.txt/;
#    $f_lib = "02";
    system "awk \'{\$1=\"$f_lib\_\"\$1; print}\'  $f  > tmp/tmp_$f_lib.txt";
}

system "cat tmp/tmp*.txt | sort -k 4 -n |sed -e 's/\\s/\\t/g' > tmp_input.txt";

my %Read = ();
my %Begin = ();
my $count = 1;
my $max_End = 1;
open F, "tmp_input.txt" or die "Cannot open file: tmp_input.txt, $!";
while(<F>){
    chomp;
    my ($t_ID, $t_Begin, $t_End, $t_Strand) = (split /\s+/)[0, 3, 4, 5];
    my $t_Length = $t_End - $t_Begin + 1;
    $max_End = ($max_End > $t_End)?$max_End:$t_End;
    my $count_tag = sprintf"C%05d",$count;
    $Read{$count_tag} = $_;
    $Begin{$t_Begin}->{$count_tag} = $_;
    $count ++;
}
close F;

################################################################################
# Merge all tags
#open OUT, ">merge.txt" or die;
my %del = ();
my $tag = 1;
foreach my $in (sort keys %Read){
    next if exists $del{$in};
    my ($in_ID, $in_Begin, $in_End, $in_Strand) = (split /\s+/, $Read{$in})[0, 3, 4, 5];
    my ($Range_L, $Range_R) = &RangeGap($in_Begin, 1000);
    my ($new_Begin, $new_End) = ($in_Begin, $in_End);
    my @tank = ();
    for(my $i=$Range_L; $i<=$Range_R; $i++){
        next unless exists $Begin{$i};
        foreach my $b (sort keys %{$Begin{$i}}){
            my ($p_ID, $p_Begin, $p_End, $p_Strand) = (split /\s+/, $Read{$b})[0, 3, 4, 5];
# Strand            
            next unless ($in_Strand eq $p_Strand);
# Overlap
            next unless ($new_Begin <= $p_End && $new_End >= $p_Begin);
# >>>            
            my $overlap = &CalOverlap($new_Begin, $new_End, $p_Begin, $p_End);
            next unless ($overlap > 1);
            my $overlap_p1 = $overlap/($p_End-$p_Begin+1);
            my $overlap_p2 = $overlap/($new_End-$new_Begin+1);
            next unless ($overlap >= 10 || $overlap_p1 >= 0.5 || $overlap_p2 >= 0.5);
# <<<
            ($new_Begin, $new_End) = &merge2to1($new_Begin, $new_End, $p_Begin, $p_End);
            $Range_R = $new_End + 1000;
            push @tank, $Read{$b};
            $del{$b} = 1;
        }
    }

    my $new_Length = $new_End - $new_Begin + 1;
    my $new_ID = sprintf"Seed%04d", $tag;
    $tag ++;
    my $new_Info = join"\t",("$new_ID", "exp", $new_Length, $new_Begin, $new_End, $in_Strand);
    print $new_Info,"\n";
#   print join"\n", (">$new_Info", @tank),"\n";
}
#close OUT;

################################################################################
sub RangeGap{
    my ($begin, $length) = @_;
    my $Left = $begin - $length;
    my $Right = $begin + $length;
    $Left = ($Left < 1)?1:$Left;
    $Right = ($Right > $max_End)?$max_End:$Right;
    return($Left, $Right);
}

sub merge2to1{
    my ($a_begin, $a_end, $b_begin, $b_end) = @_;
    my $new_begin = ($a_begin < $b_begin)?$a_begin:$b_begin;
    my $new_end   = ($a_end   < $b_end  )?$b_end:$a_end;
    return($new_begin, $new_end);
}

sub CalOverlap{
    my ($c_begin, $c_end, $d_begin, $d_end) = @_;
    my $c_length = $c_end - $c_begin + 1;
    my $d_length = $d_end - $d_begin + 1;
    my $overlap = 1;
    if($d_begin < $c_begin){
        $overlap = $d_end - $c_begin + 1;
        $overlap = $c_length if ($overlap > $c_length);
    }elsif($d_begin < $c_end){
        $overlap = $c_end - $d_begin + 1;
        $overlap = $d_length if ($overlap > $d_length);
    }else{
        $overlap = -1;
    }
    return $overlap;
}
