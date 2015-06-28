#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;

sub help{
    print STDERR <<EOF
Usage: perl getRscript_v2.pl  -i input.txt  -t 1 -m 104 -n 30
Options:
    -i <file>   : Input file. eg: match_mRNA.txt.
    -t          : Mapping types, [1] only PosCov, [3] PosCov,3/5 end.
    -m          : Totle lines in match_mRNA.txt. <Position>
    -n          : Total lines in match_mRNA.list. <sRNA>
    -h          : Show this help.

Note:
1. "104" indicated the total lines in input.txt
2. "30" indicated the totoal lines in input.match 

Description:
Mapping sRNAs that have overlap with input sequence.
EOF
}

my ($infile, $type, $max_line, $max_sRNA, $help);
&GetOptions(
            "i=s"   =>  \$infile,
            "t=s"   =>  \$type,
            "m=s"   =>  \$max_line,
            "n=s"	=>	\$max_sRNA,
            "h"     =>  \$help,
           );

if($help or not $infile or not $max_line or not $max_sRNA or not $type){
	&help();
	exit(1);	
}

unless($type == 1 || $type == 3){
    &help();
    print "Change -t to 1 or 3.\n";
    exit(1);
}

######################################
# 1. Create title and directory files.
######################################

my @types = ("PosCov", "5end", "3end");
my (@titles, @dirs);

open IN,"<$infile" or die;
while(<IN>){
    chomp;
    my @ps  = split(/\s+/);
    my $id  = $ps[0];
    my $lab = 1;

    foreach my $t (@types){
        next if($lab > $type);
        push @titles,("$id\_$t");
        push @dirs,("Cov\_$t\/$id\_$t\.txt");
        $lab ++;
    }
}

close IN;

my $title = join"\n", @titles;
my $dir   = join"\n", @dirs;

open O1,"> title.txt" or die;
print O1 $title,"\n";
close O1;

open O2,"> dir.txt" or die;
print O2 $dir,"\n";
close O2;

######################################
# 2. Create R script
######################################
my $R_pdf = $infile;
$R_pdf =~ s/\.txt/\.pdf/g;

my $R_line1="pdf\(\"$R_pdf\"\)\n".'layout(matrix(1:12, 3, 2, byrow=TRUE))';
my $R_line2='
sRNAs  <-  read.table(file="title.txt")
file  <-  read.table(file="dir.txt")
sRNAs <-  as.character(sRNAs$V1)
file  <-  as.character(file$V1)
';

my $k = ($type == 1)?"k = i":"k = ceiling(i/3)";

my $R_line3= "pos  \<\- scan\(file\=\"$infile\"\,  what\=as\.list\(rep\(\"\"\,$max_line\)\)\)".'
for(i in 1:length(sRNAs)){
	data <- read.table(file[i])
	ymax <- ceiling(max(data$V2, data$V3, data$V4, data$V5)/1000)*1000
	plot(seq(min(data$V1),max(data$V1), length=3), seq(1, ymax, length=3), 
		type="n", xlab="", ylab="", xaxt="n", yaxt="n")
# The reads coverage in four libraries.
	lines(data$V1, data$V2, pch=20, cex=2, col=2)
	lines(data$V1, data$V3, pch=20, cex=2, col=3)	
	lines(data$V1, data$V4, pch=20, cex=2, col=4)
	lines(data$V1, data$V5, pch=20, cex=2, col=5)

# For PosCov, 5, 3 end  or only PosCov one line.
'.$k.'
#    k=i
# yaxis.
    ymid=ymax/2;       y2=(ymax + ymid)/2 
# The genome coordination of annotated sRNA/mRNA
    b_m = as.numeric(pos[[10]][k]);    e_m = as.numeric(pos[[11]][k]);    strand=as.character(pos[[12]][k])
    if(strand == "+")  tag  = 2  else  tag  = 1
# Arrow-1 to 30, for the sRNA identified in this work.
'."    for\(m in 1\:$max_sRNA\)\{".'
        begin = as.numeric(pos[[(m+3)*3+1]][k]);   end = as.numeric(pos[[(m+3)*3+2]][k])
        if((begin == b_m) && (end == e_m)) yx <- ymid else yx <- y2     # check whether this is the annotated sRNA/mRNA.
        arrows(begin, yx, end, yx, angle=30, length=0.05, cod=tag,
                       col="black", lwd=1)
        lines(c(begin, begin), c(0,yx), col="grey70", lwd=0.5, lty=2)
        lines(c(end,   end),   c(0,yx), col="grey70", lwd=0.5, lty=2) 
       }

# Arrow-1, the annotated sRNA/mRNA. 
    arrows(b_m, ymid,e_m, ymid, angle=30, length=0.1, code=tag,
            col="black", lwd=2)
    lines(c(b_m, b_m), c(0,ymid), col="grey70", lwd=0.3, lty=2)
    lines(c(e_m, e_m), c(0,ymid), col="grey70", lwd=0.3, lty=2)

# Draw the axis
	axis(side=1, seq(min(data$V1), max(data$V1), length=3), tcl=-0.2, labels=FALSE)
	axis(side=2, seq(0,ymax, length=3), tcl=-0.2, labels=FALSE)
	mtext(seq(min(data$V1), max(data$V1), length=3), side=1, las=1, line=0.5,
		at=seq(min(data$V1), max(data$V1), length=3), 
		cex=0.7, font.axis=1)
	mtext(seq(0, ymax/1000, length=3), side=2, las=1, line=0.5,
		at=seq(0,ymax, length=3), 
		cex=0.7, font.axis=1)
# Draw the axis labels
	mtext("Number of reads /X1000", side=2, las=0, line=1.8, cex=0.5, font=1)
	mtext("Genome coordinate", side=1, las=1, line=1.8, cex=0.5, font=1)
# The legend	
	legend("topright", c("18-40 nt", "40-80 nt", "80-140 nt", ">140 nt"),
		col=c(2,3,4,5), text.col = "black", cex=0.7,
		pch = c(20,20,20,20))		
# The title
    Name <- sub("_PosCov","", sRNAs[i])
	title(Name, font.main=2, adj=0.5, cex.main=1.3, line=0.5)
}
dev.off()
';

my $R_out = $infile;
$R_out =~ s/\.txt/2pdf\.R/;
open OUT,">$R_out" or die;
print OUT "$R_line1\n$R_line2\n$R_line3\n";
close OUT;
