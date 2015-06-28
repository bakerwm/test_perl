#!/usr/bin/perl -w
use strict;

sub help{
    print STDERR <<EOF
Usage: perl  getRscript.pl  infile  1

Note:
1. Format of infile should be: 
ID* Exp Length  Begin*  End*    Strand*  (* required)

2. The second parameter should be: 1 or 3.
"1" for PosCov statement only
"3" for PosCov, 5end and 3end statement.

EOF
}

my $infile = shift or die &help;
my $num    = shift or die &help;

unless($num == 1 || $num == 3){
    &help();
    print "Change the second parameter to 1 or 3.\n";
    exit(1);
}

######################################
# 1. Create title and directory files.
######################################

my @types = ("PosCov", "5end", "3end");
my (@titles, @dirs);

open F,"< $infile" or die;
while(<F>){
    chomp;
    my @ps  = split(/\t/);
    my $id  = $ps[0];
    my $lab = 1;

    foreach my $t (@types){
        next if($lab > $num);
        push @titles,("$id\_$t");
        push @dirs,("Cov\_$t\/$id\_$t\.txt");
        $lab ++;
    }
}
close F;

my $title = join"\n", @titles;
my $dir   = join"\n", @dirs;

open O1,"> title.txt" or die;
print O1 $title,"\n";
close O1;

open O2,"> dir.txt" or die;
print O2 $dir,"\n";
close O2;

#######################################
# 2. Generate R script.
#######################################
my $R_pdf =  $infile;
   $R_pdf =~ s/\.txt/\.pdf/g;

my $R_line1="pdf\(\"$R_pdf\"\)\n".'layout(matrix(1:6, 3, 2))';

my $R_line2='
sRNAs <-  read.table(file="title.txt")
file  <-  read.table(file="dir.txt")
sRNAs <-  as.character(sRNAs$V1)
file  <-  as.character(file$V1)
';

my $k = ($num == 1)?"k = i":"k = ceiling(i/3)";
my $R_line3='pos  <- read.table(file = "'.$infile.'")
for(i in 1:length(sRNAs)){
	data <- read.table(file[i])
	ymax <- ceiling(max(data$V2, data$V3, data$V4, data$V5)/1000)*1000
	plot(seq(min(data$V1),max(data$V1), length=3), seq(1, ymax, length=3), 
		type="n", xlab="", ylab="", xaxt="n", yaxt="n")

	lines(data$V1, data$V2, pch=20, cex=2, col=2)
	lines(data$V1, data$V3, pch=20, cex=2, col=3)	
	lines(data$V1, data$V4, pch=20, cex=2, col=4)
	lines(data$V1, data$V5, pch=20, cex=2, col=5)
	
	ymid <- ymax/2
'."\t$k".'
#   k=ceiling(i/3);	
	if(pos$V6[k] == "+") tag <- 2 else tag <- 1
# arrow 
    arrows(pos$V4[k],ymid,pos$V5[k],ymid, angle=30, length=0.1, code=tag,
		col="black", lwd=2)	
	lines(c(pos$V4[k],pos$V4[k]),c(0,ymid), col="black", lwd=0.5, lty=2)
	lines(c(pos$V5[k],pos$V5[k]),c(0,ymid), col="black", lwd=0.5, lty=2)		

# The legend.	
	axis(side=1, seq(min(data$V1), max(data$V1), length=3), tcl=-0.2, labels=FALSE)
	axis(side=2, seq(0,ymax, length=3), tcl=-0.2, labels=FALSE)
	
	mtext(seq(min(data$V1), max(data$V1), length=3), side=1, las=1, line=0.5,
		at=seq(min(data$V1), max(data$V1), length=3), 
		cex=0.7, font.axis=1)
	mtext(seq(0, ymax/1000, length=3), side=2, las=1, line=0.5,
		at=seq(0,ymax, length=3), 
		cex=0.7, font.axis=1)
	mtext("Number of reads /X1000", side=2, las=0, line=1.8, cex=0.5, font=1)
	mtext("Genome coordinate", side=1, las=1, line=1.8, cex=0.5, font=1)
	
	legend("topright", c("18-40 nt", "40-80 nt", "80-140 nt", ">140 nt"),
		col=c(2,3,4,5), text.col = "black", cex=0.7,
		pch = c(20,20,20,20))		
# The title
    Name <- sub("_PosCov","", sRNAs[i])
    title(Name, font.main=2, adj=0.5, cex.main=1.3, line=0.5)
}
dev.off()';

my $R_out =  $infile;
   $R_out =~ s/\.txt/2pdf\.R/;

open  OUT, ">$R_out" or die;
print OUT  "$R_line1\n$R_line2\n$R_line3\n";
close OUT;
