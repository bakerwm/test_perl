#!/usr/bin/perl
#use File::Basename;
#use File::Path;
#use Cwd;
#use lib "/share/raid1/genome/biosoft/ViennaRNA-1.7.1/Perl";
#use lib "/share/raid1/genome/biosoft/ViennaRNA-1.7.1/Perl/blib/arch";
#use lib "/share/raid1/genome/biosoft/ViennaRNA-1.7.1/Perl/blib/lib";
#use RNA; # Vienna RNA package perl Interface

use lib "/home/wangming/Documents/ViennaRNA-2.0.7/Perl";
use lib "/home/wangming/Documents/ViennaRNA-2.0.7/Perl/blib/arch";
use lib "/home/wangming/Documents/ViennaRNA-2.0.7/Perl/blib/lib";
use RNA;


$file=shift;
unless(-e $file)
{
	print "perl creat_rnastructure.pl mireap-*.table\n";
	exit;
}
open IN,$file or die $!;
while(<IN>)
{
	chomp;
	next if($.==1);
	@d=split(/\t/,$_);
	($string, $structure, $ffname)=($d[8],$d[9],$d[0]);
	RNA::PS_rna_plot($string, $structure, $ffname);
	system("convert $ffname $ffname\.png");system("rm $ffname");
}
close IN;
exit;
