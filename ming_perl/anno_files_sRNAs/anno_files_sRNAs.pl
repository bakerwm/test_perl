#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd qw(cwd abs_path);
use File::Path qw(make_path remove_tree);
use File::Basename qw(basename dirname);
use File::Spec::Functions qw(catfile catdir);
use Getopt::Long;

###########################################################################
# anno_file_sRNAs.pl                                                      #
#                                                                         #
# This script is designed to add annotation to query files (db) according #
# to the comparion with sRNA files by genome coordinates. (BED format)    #
#                                                                         #
###########################################################################

sub usage {
    die("
Usage: $0 [options] <db.txt>

Options: -i  <STR>          : The directory that contains sRNA.txt files,
                              or a list of directories, each line contains only one directory.
                              The directory could be: /bam_tags/*_sRNA.txt
                                                  or: /bam_tags/scan_cutoff/update/output/*_sRNA.txt
         -t  <STR>          : The type of tags: (sRNA|UTR|IM) to compare with db file.
         -n  <STR>          : The name of the chromosome. default: [chr]
         -o  <STR>          : The output folder, save the results. [comp_db]
         -h                 : Show this help

Example:
perl $0 -i bam_tags -o comp_db known_sRNA.txt
\n");
}

##############################
# External scripts/tools     #
##############################
my $tmp_sort2bed = "/home/wangming/work/bin/tmp.sort2bed.pl";
my $update       = "/home/wangming/work/bin/temp/update_table.pl";
my $bedtools     = "/home/wangming/Documents/bedtools2-2.23.0/bin/bedtools";

###
anno_file();

###
sub anno_file {
    my $out_dir  = abs_path('comp_db');
    my $type     = 'sRNA';
    my $chr_name = 'chr';
    my $in_dir;
    my $help;
    GetOptions(
        'i|input=s'     => \$in_dir,
        'o|output=s'    => \$out_dir,
        't|type=s'      => \$type,
        'n|name=s'      => \$chr_name,
        'h|help'        => \$help,
    );

    ###
    usage() if($help);
    usage() if(@ARGV == 0 && -t STDIN);
    die("[-i] reuqired, see -h\n") if( ! defined $in_dir );
    die("[-t $type] unknown type, see -h\n") if( ! $type =~ /^(sRNA|UTR|IM)$/i);
    my $db_file  = abs_path($ARGV[0]);

    ### Prepre working dirs
    my ($input_dir, $db_dir, $comp_dir, $output_dir) = prep_wd($out_dir);

    ### Parsing the tag_dirs (single dir, or a file contain dirs (sorted))
    my @tag_files = @{ parse_tags($in_dir, $type, $input_dir) }; # output BED files

    ### Convert the db files to BED format
    my $db_bed = convert_sort2bed($db_file, $db_dir, $chr_name);
    my ($db_bed_name) = basename($db_bed) =~ /(^.*)\.bed/;
    ### compare db and tag_files
    my @comp_outs = ();
    for my $t (sort @tag_files) {
        my ($t_name) = basename($t) =~ /(^.*)\.txt/;
        my $t_bed    = convert_sort2bed($t, $input_dir, $chr_name);
        my $comp_bedname = $db_bed_name . "_vs_" . $t_name . ".bed";
        my $comp_bed     = catfile($comp_dir, $comp_bedname);
        my $cmd = "$bedtools intersect -s -wo -a $db_bed -b $t_bed > $comp_bed";
        if( system("$cmd") == 0) {
            push @comp_outs, $comp_bed;
        }else {
            die("[Failed to compare files] \n$t \nand \n$db_bed\n");
        }
    }
    
    ### Add annotation to db_file
    my @cmd_paras = ();
    for my $p (sort @comp_outs) {
        $p = abs_path($p);
        push @cmd_paras,  "perl $update -add -n 4 -A 10 -B 10 -d " . $p;
    }
    my $cmd_line  = join(' | ', "cat $db_file", @cmd_paras);
    my $anno_file = catfile($output_dir, basename($db_file));
    $anno_file    =~ s/\.txt/.anno.$type.txt/;
    if(system("$cmd_line \> $anno_file") == 0){
    }else {
        die("Fail to add annotation\n");
    }

    ### Statistics anno_file
    my $anno_stat = $anno_file;
    $anno_stat    =~ s/\.txt/.stat/;
    stat_annotation($anno_file, $anno_stat, $db_file);
    print "Finish annotation.\n";
}

### Statistic annotation file
sub stat_annotation {
    my $in  = $_[0];
    my $out = $_[1];
    my $db  = $_[2];
    ### parse the width of db file
    my $tab_sum  = 0;
    my $line_num = 0;
    open my $fh_db, "< $db" or die "Cannot open $db, $!\n";
    while(<$fh_db>) {
        chomp;
        next if(/^\#|^\s*$/);
        my @tabs = split/\s+/;
        $tab_sum += @tabs;
        $line_num ++;
    }
    close $fh_db;
    my $db_width = sprintf"%d", $tab_sum / $line_num;
    my $db_range_r = $db_width - 1;
    ### check input file
    open my $fh_in, "< $in" or die "Cannot open $in, $!\n";
    open my $fh_out, "> $out" or die "Cannot open $out, $!\n";
    while(<$fh_in>) {
        chomp;
        my @tabs = split/\s+/;
        my $range_r = @tabs - 1;
        my @stat = @tabs[0..$db_range_r];
        for my $t (@tabs[$db_width..$range_r]) {
            if($t eq '-') {
                push @stat, '0';
            }elsif($t =~ /\w+/) {
                push @stat, '1';
            }else {
            }
        }
        print $fh_out join("\t", @stat) . "\n";
    }
    close $fh_in;
    close $fh_out;
}

### Convert the SORT (TAB) file to BED format
sub convert_sort2bed {
    my $in_txt   = $_[0];
    my $bed_dir  = $_[1];
    my $chr_name = $_[2];
    ###
    my ($bed_filename) = basename($in_txt) =~ /(^.*)\.txt$/;
    my $bed_file = catfile(abs_path($bed_dir), $bed_filename . '.bed');
    my $cmd = "perl $tmp_sort2bed -a sort -b bed -n $chr_name -x 0 $in_txt > $bed_file";
    if( system($cmd) == 0 ) {
    }else {
        print "Failed to execute [sort2bed] for: $in_txt, $!\n";
    }
    return($bed_file);
}

### Parse specific files in txt format (sRNA|UTR|IM)
sub parse_tags {
    my $in      = $_[0];
    my $type    = $_[1];
    my $out_dir = $_[2];
    my @hit_dirs = ();
    if( -f $in) {
        ### parse the dirs in the file
        open my $fh_in, "< $in" or die "Cannot open $in, $!\n";
        while(<$fh_in>) {
            chomp;
            next if(/^\#|^\s*$/);
            my $tab = (split/\s/,$_)[0];
            push @hit_dirs, $tab;
        }
        close $fh_in;
        die("[Failed to parse dirs] $in \n") if(@hit_dirs < 1);
    }elsif( -d $in ) {
        push @hit_dirs, $in;
    }else {
        die("[Need a file or directory] $in\n");
    }
    ### parse the txt files
    my @hit_files = ();
    for my $d (sort @hit_dirs) {
        ### search specific files
        my $alt_dir   = catdir(catdir(catdir($d, 'scan_cutoff'), 'update'), 'output');
        my $check_f1  = 0;
        my $check_f2  = 0;
        my @tag_files1 = glob"$d\/*$type\.txt";
        my @tag_files2 = glob"$alt_dir\/*$type\.txt";
        my @tag_files  = ();
        if(@tag_files1 >= 1) {
            $check_f1 ++;
        }elsif(@tag_files2 >= 1) {
            $check_f2 ++;
        }else {
            die"Failed to find [$type] files in $d";
        }
        ###
        if($check_f1 >= 1) {
            @tag_files = @tag_files1;
        }elsif($check_f2 >= 1) {
            @tag_files = @tag_files2;
        }else {
            #
        }
        ### Prepare input files and 
        my $tag_file  = $tag_files[0];
        my $d_name    = basename($d);
        my $tag_newname = $d_name . "_" . basename($tag_file);
        my $tag_new     = catfile($out_dir, $tag_newname);
        link $tag_file, $tag_new;
#        my $tag_bed     = covert_sort2bed($tag_new, $out_dir, $chr_name);
#        push @hit_files, $tag_bed;
        push @hit_files, $tag_new;
    }
    return(\@hit_files);
}


### Prepare working directories
sub prep_wd {
    my $out_dir = $_[0];
    make_path($out_dir) if( ! -d $out_dir );
    my $sub_input  = catdir($out_dir, 'input');
    my $sub_db     = catdir($out_dir, 'db');
    my $sub_comp   = catdir($out_dir, 'comp');
    my $sub_output = catdir($out_dir, 'output');
    ###
    make_path($sub_input) if( ! -d $sub_input );
    make_path($sub_db) if( ! -d $sub_db);
    make_path($sub_comp) if( ! -d $sub_comp);
    make_path($sub_output) if( ! -d $sub_output);
    ###
    return($sub_input, $sub_db, $sub_comp, $sub_output);
}

### END ###
