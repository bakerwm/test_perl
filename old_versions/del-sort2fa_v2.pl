### !--------------------------------------------------
### replaced by: sort2bed.pl

#!/usr/bin/perl -w
use warnings;
use strict;
use File::Spec::Functions qw(catdir catfile);
use File::Basename qw(basename dirname);
use Getopt::Long;
use Pod::Usage;

use Data::Dumper;

my $infile = '';
my $db = '';
my $refSeq = '';
my $outfile = '';
my $help = '';
GetOptions('i|input=s' => \$infile,
           'd|database=s' => \$db,
           'f|reference=s' => \$refSeq,
           'o|outfile=s' => \$outfile,
           'h|help' => \$help
           ) or pod2usage(-verbose => 0);

pod2usage(-verbose => 1) if($help);
pod2usage(-message => "$0: [-i|--input $infile] Not found.") unless -e $infile;

## Read infile
my %Strain = ();
my %Info = ();

open F, $infile or die "$!";
while(<F>){
    chomp;
    my @tabs = split /\t/;
    warn "Line$. do not have six columns\n" if(@tabs < 6);
    my ($id, $strain) = @tabs[0,1];
    $Strain{$strain} = 1;
    $Info{$id} = $_;
}
close F;

# Retrieve genomes from $db
my $ref_dir = '/share/wangming/reference_genomes/MTB_complex';

my $outfa = $infile;
$outfa =~ s/\.txt/.fa/;
$outfile = $outfa unless($outfile);
open OUT, ">$outfile" or die "$!";
my @outSeqs = ();

if($refSeq) {
    pod2usage(-message => "File not found: $refSeq") unless -e $refSeq;
    my $genome = &ReadFA($refSeq);
    foreach my $i (keys %Info) {
        my ($id, $Start, $End, $Strand) = (split /\t/, $Info{$i})[0,3,4,5];
        my $i_Seq = &txt2fa($genome, $Start, $End);
        $i_Seq = &RevComp($i_Seq) if($Strand eq '-');
        push @outSeqs, (">$id", $i_Seq);
    }
    my $out = join "\n", @outSeqs;
    print OUT $out, "\n";
} else{
    pod2usage(-message => "Need input the dir of reference: -d") unless ($db);
# Prepare the reference genome
    my %gm = ();
    my $count = 1;
    foreach my $s (keys %Strain){
        my $s_fna = catfile($db, "$s\.fna");
        my $s_fa = catfile($db, "$s\.fa");
        my $s_file = '';
        $s_file = $s_fna if(-e $s_fna);
        $s_file = $s_fa if(-e $s_fa);
        pod2usage(-message => "Ref file not found: [$s], at [$db]") unless -e $s_file;
        $gm{$s} = &ReadFA($s_file);
        $count ++;
    }
    foreach my $i (keys %Info) {
        my ($id, $Strain, $Start, $End, $Strand) = (split /\t/, $Info{$i})[0,1,3,4,5];
        my $genome = $gm{$Strain};
        my $i_Seq = &txt2fa($genome, $Start, $End);
        $i_Seq = &RevComp($i_Seq) if($Strand eq '-');
        push @outSeqs, (">$id", $i_Seq);
    }
    my $out = join "\n", @outSeqs;
    print OUT $out, "\n";
}
close OUT;

# Subroutines 
sub ReadFA{
    my $fa = shift(@_);
    open F, $fa or die "$fa, $!";
    $/ = "\>";
    my $faSeq = '';
    while(<F>){
        s/\>//;
        my @lines = split /\n/;
        next if @lines < 2;
        my $id = shift(@lines);
        $faSeq = join '', @lines; 
    }
    close F;
    $/ = "\n";
    return $faSeq;
}

sub txt2fa{
    my ($genome, $start, $end) = @_;
    my $length = $end - $start + 1;
    $start -= 1;
    my $fa = substr($genome, $start, $length);
    return $fa;
}

sub RevComp{
    my $infa = shift(@_);
    $infa =~ tr/ATCGatcg/TAGCtagc/;
    $infa = reverse($infa);
    return $infa;
}

__END__

=head1 NAME

<sort2fa.pl> -- transform the txt to fasta file.

=head1 SYNOPSIS

# From single ref

perl sort2fa.pl -i input.txt -f ref.fa -o output.fa

# From multiple ref

perl sort2fa.pl -i input.txt -d database -o output.fa

=head1 OPTIONS

=over 8

=item B<-i> B<--input>

Input the sort file, with at least 6-col, * are required.
<ID*> <Strain> <Length> <Start*> <End*> <Strand*>

=item B<-d> B<--db>

The dir that contain the referenece files.

=item B<-f> B<--reference>

For the sort text file from the same reference (Single reference).

=item B<-o> B<--outfile>

Write the result fasta file to the file. default [input.fa]

=head1 B<-h>

Print this help.

=back

=cut
