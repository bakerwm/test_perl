#!/usr/bin/perl
#!/bin/sh
use warnings;
use strict;
use Getopt::Long;

my ($help,$input,$output,$fdr,$log2,$nameA,$nameB,$Head);
GetOptions(
	"help"=>\$help,	
	"input=s"=>\$input,
	"output=s"=>\$output,
	"fdr=f"=>\$fdr,
	"log2=f"=>\$log2,
	"nameA=s"=>\$nameA,
	"nameB=s"=>\$nameB,
	"Head"=>\$Head,
);
if (!defined $input || !defined $nameA || !defined $nameB || defined $help) {
	&usage();
	exit 1;
}

$output ||=$input.'.pdf';
$fdr ||=0.001;
$log2 ||=1;

my $tmpxls=$input.'.tmp';
open IN,$input or die "can't open the input file $input";
open OUT,'>',$tmpxls or die "can't open the out tmp xls file $tmpxls";
if (defined $Head) {
	<IN>;
}
while (my $line=<IN>) {
	chomp $line;
	my @tab=split /\t/,$line;
	if ($tab[4]==0) {
		$tab[4]=0.001;
	}
	if ($tab[5]==0) {
		$tab[5]=0.001;
	}
	my $out= join "\t",@tab;
	print OUT $out."\n";
}
close IN;
close OUT;

my $Rsh=<<RSH;
pdf(file ="$output");
data <- read.table ("$tmpxls")
plot (log10(data[abs(data[,7])<$log2 | data[,10]>$fdr,5]),log10(data[abs(data[,7])<$log2 | data[,10]>$fdr,6]), xlab = "log10($nameA RPKM)", ylab = "log10($nameB  RPKM)", col = "blue", main = "Expression Level $nameA vs $nameB", pch = ".", cex = 2)
x <- data[abs(data[,7]) >= $log2 & data[,10] <= $fdr ,5]
y <- data[abs(data[,7]) >= $log2 & data[,10] <= $fdr ,6]
for (i in 1 : length(x)) {
	if (y[i] > x[i])
		points(log10(x[i]), log10(y[i]), col = "red", pch = ".", cex = 2)
	else
		points(log10(x[i]), log10(y[i]), col = "green", pch = ".", cex = 2)
}
legend("topleft", c("up-regulated", "down-regulated", "Not DEGs"), fill = c("red", "green", "blue"), xjust = 0, title = paste("FDR <= ", $fdr, " AND |log2Ratio| >= ", $log2))
save.image()
dev.off()
RSH

my $Rfile=$output.'.R';
open RDRAW,'>',$Rfile or die "can't open the R file for $output ";
print RDRAW $Rsh;
close RDRAW;

system "/share/raid1/genome/bin/Rscript  $Rfile";
my $outpng=$output;
$outpng=~s/\.pdf$/\.png/;
system "/usr/bin/convert -density 98 $output $outpng ";

sub usage{
	print << "USAGE";
name
	plotDiffExp.pl
descripyion
	draw figure of diff expression
	-input  (str) the input file \$5,A_RPKM ; \$6,B_RPKM ; \$7 log2(\$6/\$5) ; \$10 FDR ;
	-output (str) the pdf of output 
	-fdr    (f)   default= 0.001
	-log2   (f)   default= 1
	-nameA  (str) species A  name  
	-nameB  (str) species B  name  
	-Head         the input with head ( 1 line)
	-help         help
author
	nixiaoming   nixiaoming\@genomics.cn
version
	1.0    2010-1-5 15:40
example	
	perl plotDiffExp.pl -h
	perl plotDiffExp.pl -input liver-VS-spleen.GeneDiffExp.xls -Head -output liver-VS-spleen.pdf -fdr 0.001 -log2 1 -nameA liver -nameB spleen 

USAGE
}
