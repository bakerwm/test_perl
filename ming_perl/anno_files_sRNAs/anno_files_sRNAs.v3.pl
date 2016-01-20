#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd qw(cwd abs_path);
use File::Path qw(make_path remove_tree);
use File::Basename qw(basename dirname);
use File::Spec::Functions qw(catfile catdir);
use Getopt::Long;
use Data::Dumper;

###########################################################################
# anno_file_sRNAs.pl                                                      #
#                                                                         #
# This script is designed to add annotation to query_files (db) according #
# to the comparion with sRNA files by genome coordinates. (BED format)    #
#                                                                         #
# 1. Add the ids of query_files to db file by genome coordinates. (*.txt) #
# 2. Calculate the distance between tags in db_file and queryfile in both #
#    5' and 3' ends. (*53ends.txt)                                        #
# 3. Check whether the tags in db were found in each query condition,     #
#    only report 1 or 0, represent yes or no.                             #
#                                                                         #
# 2015-12-28 Wang Ming wangmcas(AT)gmail.com                              #
###########################################################################

##############################
# External scripts/tools     #
##############################
my $tmp_sort2bed = "/home/wangming/work/bin/tmp.sort2bed.pl";
my $update       = "/home/wangming/work/bin/temp/update_table.pl";
my $bedtools     = "/home/wangming/Documents/bedtools2-2.23.0/bin/bedtools";

###
anno_file();

##############################
#   Subroutines              #
##############################
sub anno_file {
    my $out_dir  = abs_path('comp_db');
    my $type     = 'sRNA';
    my $chr_name = 'chr';
    my $db_file;
    my $help;
    GetOptions(
        'i|input=s'     => \$db_file,
        'o|output=s'    => \$out_dir,
        't|type=s'      => \$type,
        'n|name=s'      => \$chr_name,
        'h|help'        => \$help,
    );
    ###
    usage() if($help);
    usage_short() if(@ARGV == 0 && -t STDIN);
    die("[-i] reuqired, see -h\n") if( ! defined $db_file );
    die("[-t $type] unknown type, see -h\n") if( ! $type =~ /^(sRNA|UTR|IM)$/i);
    my $in_dir  = abs_path($ARGV[0]);
    ### Prepre working dirs
    my ($input_dir, $db_dir, $comp_dir, $output_dir) = prep_wd($out_dir);
    ### Parsing the tag_dirs (single dir, or a file contain dirs (sorted))
    my @tag_files = @{ parse_tags($in_dir, $type, $input_dir) }; # output BED files
    my $anno_lib  = catfile($output_dir, $type . '_anno_list.lib'); ### save the list of selected tag_files
    save_tag_lists($anno_lib, \@tag_files);
    ### Convert the db files to BED format
    my $db_bed = convert_sort2bed($db_file, $db_dir, $chr_name);
    my ($db_bed_name) = basename($db_bed) =~ /(^.*)\.bed/;
    ### Compare db and tag_files
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
    ### Evaluation the 5'/3' distance
    @tag_files = sort(@tag_files);
    my $check_53ends_file = $anno_file;
    $check_53ends_file    =~ s/\.txt$/.53ends.txt/;
    stat_53ends($anno_file, \@tag_files, $check_53ends_file, $db_file);
    print "Finish annotation.\n";
}

### Prepare work directories
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
            die"Failed to find [$type] files in $d\n";
        }
        ###
        if($check_f1 >= 1) {
            @tag_files = @tag_files1;
        }elsif($check_f2 >= 1) {
            @tag_files = @tag_files2;
        }else {
            #next;
        }
        ### Prepare input files and 
        my $tag_file  = $tag_files[0];
        my $d_name    = basename($d);
        my $tag_newname = $d_name . "_" . basename($tag_file);
        my $tag_new     = catfile($out_dir, $tag_newname);
        link $tag_file, $tag_new;
        push @hit_files, $tag_new;
    }
    return(\@hit_files);
}

### Save the list of selected tag_files to lib
sub save_tag_lists {
    my $anno_lib  = $_[0];
    my @tag_files = @{$_[1]};
    open my $fh_an, "> $anno_lib" or die "Cannot open $anno_lib, $!\n";
    print $fh_an "\# Add information from the following files to db_file\n";
    my $order = 0;
    for my $t (@tag_files) {
        $order ++;
        $order = sprintf"%02d", $order;
        print $fh_an $order . "\t" . $t . "\n";
    }
    print $fh_an "\n";
    ###
    print $fh_an "\# File: *[type].txt          : The ids from the tags files" . "\n";
    print $fh_an "\# File: *[type].stat         : Convert the ids to 1 or 0 based on the ids, hit = 1, not_hit('-') = 0" . "\n";
    print $fh_an "\# File: *[tyep].53ends.txt   : Check the distance of 5'/3' of the tags to the db. each file will output 3-column: overlap:5_distance:3_distance" . "\n";
    close $fh_an;
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

### Statistic annotation file
### check the existence of tags in db_file with query files, 0 or 1
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

### check the distance between tags and db files in 5/3 ends
sub stat_53ends {
    my $anno_file = $_[0];
    my @tag_files = @{ $_[1] };
    my $check_53ends_file = $_[2];
    my $db_file           = $_[3];
    ### parse the tag_files
    my %tag_info = ();
    my $rank = 1;
    for my $b (sort @tag_files) {
        %{$tag_info{$rank}} = %{ read_file($b) };
        $rank ++;
    }
    ### check the annotation file
    my $db_width = file_width($db_file);
    open my $fh_anno, "< $anno_file" or die "Cannot open $anno_file, $!\n";
    open my $fh_53, "> $check_53ends_file" or die "Cannot open $check_53ends_file, $!\n";
    while(<$fh_anno>) {
        chomp;
        next if(/^\#|^\s*$/);
        my @tabs = split /\t+/, $_;
        my ($db_left, $db_right, $db_str) = @tabs[3..5];
        ### set the range of the array
        my $db_range_r = $db_width - 1;
        my $tab_range_r = @tabs - 1;
        my @tabs_stats = @tabs[0..$db_range_r];
        ### scan the mapping tags
        my $order = 0;
        for my $hit (@tabs[$db_width..$tab_range_r]) {
            $order ++;
            my ($overlap, $dis_5, $dis_3) = (0, '-', '-');
            if($hit eq '-') {
                push @tabs_stats, ($overlap, $dis_5, $dis_3);
            }elsif($hit =~ /\w+/) {
                my ($query_left, $query_right) = determine_ends(\%tag_info, $order, $hit);
                ($overlap, $dis_5, $dis_3) = determine_53dis($db_left, $db_right, $query_left, $query_right, $db_str);
                push @tabs_stats, ($overlap, $dis_5, $dis_3);
            }
        }
        print $fh_53 join("\t", @tabs_stats) . "\n";
    }
    close $fh_anno;
    close $fh_53;
}

### Determine the 5'/3' ends of ids
sub determine_ends {
    my %info    = %{ $_[0] };
    my $order   = $_[1];
    my $id      = $_[2];
    ###
    my @ids = (split/\,/, $id);
    my $left_end;
    my $right_end;
    my @left_ends;
    my @right_ends;
    for my $i (@ids) {
        my ($sta, $end, $str);
        if(exists $info{$order}->{$i}) {
            ($sta, $end, $str) = split/\:/, $info{$order}->{$i};
        }else {
            print "[$i] not found in annotation file.\n";
        }
        push @left_ends, $sta;
        push @right_ends, $end;    
    }
    $left_end = min_max(\@left_ends, 'min');
    $right_end = min_max(\@right_ends, 'max');
    return($left_end, $right_end);
}

### Determine the distance between query and db in both 5' and 3' ends
sub determine_53dis {
    my $db_left     = $_[0];
    my $db_right    = $_[1];
    my $query_left  = $_[2];
    my $query_right = $_[3];
    my $db_str      = $_[4];
    ###
    my $left_dis  = $db_left - $query_left;
    my $right_dis = $query_right - $db_right;
    my ($dis_5, $dis_3) = ($left_dis, $right_dis);
    if($db_str eq '-') {
        ($dis_5, $dis_3) = ($right_dis, $left_dis);
    }
    ### calculate overlap
    my @vals = ();
    push @vals, $db_right - $query_left + 1;
    push @vals, $query_right - $db_left + 1;
    push @vals, $db_right - $db_left + 1;
    push @vals, $query_right - $query_left + 1;
    my $overlap = min_max(\@vals, 'min');
    return($overlap, $dis_5, $dis_3);
}

### Cal min/max value
sub min_max {
    my @vals = @{$_[0]};
    my $type = $_[1];
    my $min = my $max = shift(@vals);
    for my $v (@vals) {
        $min = ($min < $v)?$min:$v;
        $max = ($max > $v)?$max:$v;
    }
    if($type eq 'min') {
        return($min);
    }elsif($type eq 'max') {
        return($max);
    }else {
        die("[$type] unknown type, min or max. \n");
    }
}

### Check the width of input db_file
sub file_width {
    my $in = $_[0];
    ### parse the width of db file
    my $tab_sum  = 0;
    my $line_num = 0;
    open my $fh_in, "< $in" or die "Cannot open $in, $!\n";
    while(<$fh_in>) {
        chomp;
        next if(/^\#|^\s*$/);
        my @tabs = split/\s+/;
        $tab_sum += @tabs;
        $line_num ++;
    }
    close $fh_in;
    my $db_width = sprintf"%d", $tab_sum / $line_num;
    my $db_range_r = $db_width - 1;
    return($db_width);
}

### save tag position in hash
sub read_file {
    my $in_file = $_[0];
    my %info    = ();
    open my $fh_in, "< $in_file" or die "Cannot open $in_file, $!\n";
    while(<$fh_in>) {
        chomp;
        next if(/^\#|^\s*$/);
        my ($id, $start, $end, $strand) = (split /\t/, $_)[0, 3, 4, 5];
        $info{$id} = join("\:", $start, $end, $strand);
    }
    close $fh_in;
    return(\%info);
}

sub usage_short {
    die("
Usage: perl anno_files_tags.pl -i db_file.txt  [tag_dir|tag_dir.list]

see: -h for more details
\n");

}

sub usage {
    die("
Usage: $0 [options] <tag_dir|tag_dir.list>

Options: -i  <STR>          : The db file in tab-separate format, with at least 6 fields:
                              <id*> <chr> <length> <start> <end> <strand>
         -t  <STR>          : The type of tags: [sRNA|UTR|IM] to compare with db file.
         -n  <STR>          : The name of the chromosome. default: [chr]
         -o  <STR>          : The output folder, save the results. [comp_db]
         -h                 : Show this help

Notes:   <tag_dir|tag_dir.list>
         Support the path to a <directory> and a file contains a list of <directories>
         1. tag files could be located in [dir] or [dir]/scan_cutoff/update/output/
         2. tag files have to be in the following name: *_[sRNA|UTR|IM].txt
         3. if input a list of directorie, support only one [dir] in each line

Example:
perl $0 -i bam_tags -o comp_db known_sRNA.txt
\n");
}

### END ###
