#!/usr/bin/env perl 

#################################################
# update_table.pl
# 
# manipulate two files: add/update columns
# 1. info table. contain the info ready to update
# 2. index table. contain the id index of both
#    files.
#
# Wang Ming wangmcas(at)gmail.com
# 2015-07-07
#################################################

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case); # case senstive
use Data::Dumper;

update_table();
#exit(1);

sub update_table {
    my $db    = '';
    my $cm    = 1;
    my $ca    = 1;
    my $cb    = 1;
    my $cn    = 1;
    my $cA    = 1;
    my $cB    = 1;
    my $dict  = '';
    my $ci    = 1;
    my $cj    = 2;
    my $hit_type = 0;
    my $add;
    my $replace;
    my $line;
    my $output;
    my $log;
    my $report;
    my $help;
    my $help_detail;
    GetOptions(
            'd=s'       => \$db,
            'm=i'       => \$cm,
            'a=i'       => \$ca,
            'b=i'       => \$cb,
            'n=i'       => \$cn,
            'A=i'       => \$cA,
            'B=i'       => \$cB,
            'dict|x=s'  => \$dict,
            'i=i'       => \$ci, # query id
            'j=i'       => \$cj, # db id
            'o=s'       => \$output,
            'hit|U=i'   => \$hit_type,
            'add|D'     => \$add,
            'replace|R' => \$replace,
            'line|L'    => \$line,
            'report=s'  => \$report,
            'log=s'     => \$log,
            'h'         => \$help,
            'help'      => \$help_detail,
            ) or die "see: perl $0 -h , for help\n";
    usage_full() if($help_detail);
    usage_simple() if($help);
    if(@ARGV == 0 && -t STDIN) {
        die("Usage: perl $0 [-add|-replace|-line] [-m -a -b] [-d -n -A -B] <input>

see -h or --help for more details\n");
    }
    
    # save db information table to hash
    die("[-d $db] not exists\n") if(! -e $db);
    my %df = %{read_file($db, $cn)}; # save db to hash by ids as keys

    ### ----prepare dict hash----------------------------------------------
    # convert -j to -i: ids of db to ids of query/input/STDIN
    if($dict) { # information to in.file
        my %dict_converter = %{read_file($dict, $ci, $cj)}; # i=query-id, j=db-id
        my %df_new = ();
        my $dict_status = 0;
        for my $d (keys %df) {
            next if( ! exists $dict_converter{$d} );
            $df_new{$dict_converter{$d}} = $df{$d};
            $dict_status ++;
        }
        die("ids of db (-d -n) in NOT found in dict table (-x -i)\n");
        %df = ();
        %df = %df_new;
    }
    
    ### ----set report type, change lot ...--------------------------------
    die("unknown [-hit $hit_type], expect:0, 1 or 2\n") 
        unless($hit_type == 0 || $hit_type == 1 || $hit_type == 2);
    
    ### ----START processing query/STDIN file------------------------------
    my $fh_out = *STDOUT;
    if(defined $output) {
        open $fh_out, "> $output" or die "Cannot write to $output, $!\n";
    }
    my $report_out = '';
    my $change_log = '';
    my $change_num = 0;
    if($add) {
        while(<>) {
            chomp;
            my ($new_line, $log) = add_tabs($_, $cm, $ca, \%df, $cA, $cB, $hit_type, $log);
            $change_log .= $log . "\n";
            next unless($new_line); # skip 0 
            print $fh_out $new_line . "\n";
            $change_num ++;
        }
        $report_out = "Type:\tAdd columns [$cA\-$cB] to last column STDIN/input";
    }elsif($replace) {
        while(<>) {
            chomp;
            my ($new_line, $log) = replace_fields($_, $cm, $ca, $cb, \%df, $cA, $cB, $hit_type, $log);
            $change_log .= $log . "\n";
            next unless($new_line); # skip 0
            print $fh_out $new_line . "\n";
            $change_num ++;
        }
        $report_out = "Type:\tReplace column-[$ca\-$cb] of STDIN/input by column-[$cA\-$cB] of (-d) file";
    }elsif($line) {
        while(<>) {
            chomp;
            my ($new_line, $log) = extract_lines($_, $cm, \%df, $hit_type, $log);
            $change_log .= $log . "\n";
            next unless($new_line); # skip 0
            print $fh_out $new_line . "\n";
            $change_num ++;
        }
        $report_out = "Type:\tExtract rows from [-d] file";
    }else {
        die("[-add|-line|-replace] not found, see -h\n");
        exit(1);
    }
    close $fh_out;

# redirect change log to file
    if($log) {
        open my $fh_log, "> $log" or die "Cannot write to $log, $!\n";
        print $fh_log $change_log . "\n";
        close $fh_log;
    }

# generate a report 
    if($report) {
        open my $fh_rpt, "> $report" or die "Cannot wirte to $report, $!\n";
        print $fh_rpt $report_out . "\n";
        close $fh_rpt;
    }
}

###################################
# Module-1: add (insert) tabs 
# ? support one to N ?
###################################
sub add_tabs { 
    my $in       = $_[0];
    my $m        = $_[1];
    my $a        = $_[2];
    my %df       = %{$_[3]};
    my $A        = $_[4];
    my $B        = $_[5];
    my $hit_type = $_[6];
    my $log      = $_[7];
    # 
    my $new_in = '';
    my $ch_log = '';
    my $add_width  = $B - $A + 1;
    my @add_blanks = split //, ('-' x $add_width);
    # 
    if($in =~ /^\s*$|^\#/) {
        $new_in = $in;
    }else {
        my $insert_mark = '';
        my $in_left  = $in;
        my $in_right = '';
        # split input line
        if($a > 1) { # need to add insert in line
            my @in_left_tabs  = extract_tabs($in, 1, $a);
            my @in_right_tabs = extract_tabs($in, ($a + 1), 'end');
            $in_left  = join("\t", @in_left_tabs);
            $in_right = join("\t", @in_right_tabs);
            $insert_mark = "column-$a";
        }else { # default, to the tail of line
            $in_left = $in;
            $insert_mark = 'end';
        }
        my @in_ids = extract_tabs($in, $m, $m);
        my $in_id  = shift(@in_ids);
        my @insert_tabs = '';
        ### hit type
        if($hit_type == 0) {
### support 1 => N, id mapping
            if(exists $df{$in_id}) {
#                @insert_tabs = extract_tabs($df{$in_id}, $A, $B);
                @insert_tabs = @{ support_nTabs(\@{$df{$in_id}}, $A, $B) };
            }else { # only report not hit
                @insert_tabs = @add_blanks;
            }
            $new_in = join("\t", $in_left, @insert_tabs, $in_right);
        }elsif($hit_type == 1) {
            if(exists $df{$in_id}) {
#                @insert_tabs = extract_tabs($df{$in_id}, $A, $B);
                @insert_tabs = @{support_nTabs(\@{$df{$in_id}}, $A, $B) };
                $new_in = join("\t", $in_left, @insert_tabs, $in_right);
            }else{
                $new_in = 0;
            }
        }elsif($hit_type == 2) {
            if(exists $df{$in_id}) {
                $new_in = 0;
            }else{
                @insert_tabs = @add_blanks;
                $new_in = join("\t", $in_left, @insert_tabs, $in_right);
            }
        }else {
            die("[-hit $hit_type] unknown hit type, expect 0, 1 or 2\n");
        }
        $new_in =~ s/\s+$//s; # trim the blanks in tail
        $ch_log = "Add to $add_width fields after [$insert_mark] of input file" if($log);
    }
    return($new_in, $ch_log);
}

##################################
# Module-2: replace (update) tabs
###################################
sub replace_fields {
    my $in       = $_[0];
    my $m        = $_[1];
    my $a        = $_[2];
    my $b        = $_[3];
    my %df       = %{$_[4]};
    my $A        = $_[5];
    my $B        = $_[6];
    my $hit_type = $_[7];
    my $log      = $_[8];
    ### check the width of input and output
    my $origin_width = $b - $a + 1;
    my $insert_width = $B - $A + 1;
    die("[-a -b, $a $b] and [-A -B, $A $B] differ in width, only accept same width replacement\n") 
        if($origin_width != $insert_width);
    ###
    my $new_in = '';
    my $ch_log = '';
    ### split in to 3 parts: left, mid, right
    my @in_left_tabs  = extract_tabs($in, 1, ($a -1));
    my @in_mid_tabs    = extract_tabs($in, $a, $b);
    my @in_right_tabs = extract_tabs($in, ($a + 1), 'end');
    ### get insert tabs
    if($in =~ /^\s*$|^\#/) {
        $new_in = $in;
    }else{
        my @insert_tabs = ();
        my @in_ids = extract_tabs($in, $m, $m);
        my $in_id  = shift(@in_ids);
        if($hit_type == 0) {
            if(exists $df{$in_id}) {
# support 1=>N tabs
                @insert_tabs = @{ support_nTabs(\@{$df{$in_id}}, $A, $B) };
#                @insert_tabs = extract_tabs($df{$in_id}, $A, $B);
                $ch_log = "replace [col-$a to $b] " . join("\t", @in_mid_tabs) . 
                    " by " . join("\t", @insert_tabs) if($log);
            }else {
                @insert_tabs = @in_mid_tabs;
            }
            $new_in = join("\t", @in_left_tabs, @insert_tabs, @in_right_tabs);
        }elsif($hit_type == 1) {
            if(exists $df{$in_id}) {
# support 1=>N tabs
                @insert_tabs = @{ support_nTabs(\@{$df{$in_id}}, $A, $B) };
#                @insert_tabs = extract_tabs($df{$in_id}, $A, $B);
                $ch_log = "replace [col-$a to $b] " . join("\t", @in_mid_tabs) .
                    " by " . join("\t", @insert_tabs) if($log);
                $new_in = join("\t", @in_left_tabs, @insert_tabs, @in_right_tabs);
            }else{
                $new_in = 0;
            }
        }elsif($hit_type == 2) {
            if(exists $df{$in_id}) {
                $new_in = 0;
            }else {
                @insert_tabs = @in_mid_tabs;
                $new_in = join("\t", @in_left_tabs, @insert_tabs, @in_right_tabs);
            }
        }else{
            die("[-hit $hit_type] unknown hit type, expect 0, 1 or 2\n");
        }
        $new_in =~ s/\s+$//s;
    }
    return($new_in, $ch_log);
}

###################################
# Module-3: line (extract lines)
# : default only report 
#   hit lines (hit=1)
###################################
sub extract_lines {
    my $in       = $_[0];
    my $m        = $_[1];
    my %df       = %{$_[2]};
    my $hit_type = $_[3];
    my $log      = $_[4];
    #
    my $new_in = '';
    my $ch_log = '';
    my @in_ids = extract_tabs($in, $m, $m);
    my $in_id  = shift(@in_ids);
   
    if($hit_type == 0 || $hit_type == 1) {
        if(exists $df{$in_id}){
# support 1=>N tabs
            die("multiple hits found, not allowed: $in_id\n") if(@{$df{$in_id}} > 1);            
            $new_in = shift (@{$df{$in_id}});
            $ch_log = "extract line from db" if($log);
        }else {
            $new_in = 0;
        }
    }elsif($hit_type == 2) {
        if(exists $df{$in_id}) {
            $new_in = 0;
        }else {
            $new_in = $in;
        }
    }else{
        die("[-hit $hit_type] unknown hit type, expect 0, 1 or 2\n");
    }
    $new_in =~ s/\s+$//;    
    return($new_in, $ch_log);    
}

#########################
# save file to hash (id)
# # support 1 to N.
#########################
sub read_file {
    my $in      = $_[0];
    my $col_id  = $_[1];
    my $col_val = $_[2];
    my %df = ();
    open my $fh_in, "< $in" or die "Cannot open $in, $!\n";
    while(<$fh_in>) {
        chomp;
        next if(/^\s*$|^\#/);
        my @tabs = split /\t/;
        next if($col_id > @tabs); # col-id out range of file
        if(defined $col_val) {
            if($col_val > @tabs) {
                die("[out range of file, see \%dict_converter]\n");
            }
# suuport 1=>N tabs
            push @{$df{$tabs[$col_id - 1]}}, $tabs[$col_val - 1];
#            $df{$tabs[$col_id - 1]} = $tabs[$col_val - 1];        
        }
        push @{$df{$tabs[$col_id - 1]}}, $_;
#        $df{$tabs[$col_id - 1]} = $_;
    }
    close $fh_in;
    return(\%df);
}

#########################
# Extract tabs
# : support INT, 'end'
#########################
sub extract_tabs {
    my $in    = $_[0];
    my $start = $_[1];
    my $end   = $_[2];
    ###
    my @tabs  = split /\t/, $in;
    $end = scalar(@tabs) if($end =~ /^end$/i);
    die("ERROR: need INT >= 1, check (-a, -b, -A, -B, -i, -j") if($start < 1 || $end < 1);
    die("ERROR: exceed width, check (-a, -b, -A, -B, -i, -j)\n") if($end > @tabs);
    die("ERROR: need start < end, check (-a < -b, -A < -B, -i < -j)\n") if($start > $end);
    my @out = @tabs[($start - 1)..($end - 1)];
    return @out;
}

########################
# combine each tabs by ','
# # support 1=>N tabs
########################
sub support_nTabs {
    my @lines = @{$_[0]};
    my $A     = $_[1];
    my $B     = $_[2];
    my @insert_tabs = ();
    if(@lines == 1) {
        @insert_tabs = extract_tabs($lines[0], $A, $B);    
    }elsif(@lines > 1) {
        my @temp_tabs = ();
        my $temp_width = 1;
        for my $sub_line (@lines) {
            my @sub_tabs = extract_tabs($sub_line, $A, $B);
            $temp_width  = @sub_tabs;
            push @temp_tabs, \@sub_tabs;
        }
        for(my $i = 1; $i <= $temp_width; $i ++) {
            my @i_tabs = ();
            for my $t (@temp_tabs) {
                my @t_tabs = @{$t};
                push @i_tabs, $t_tabs[$i - 1];
            }
            push @insert_tabs, join(",", @i_tabs);
        }
    }
    return(\@insert_tabs);
}


### delete?
sub guess_file_width {
# return the columns of the file
    my $in   = $_[0];
    ###
    my $width_all = 0;
    my $count  = 0;
    open my $fh_in, "< $in" or die "Cannot open $in, $!\n";
    while(<$fh_in>) {
        last if($count > 500);
        chomp;
        next if(/^\s*$|^\#/);
        my @tabs = split /\t/;
        $width_all += scalar(@tabs);
        $count ++;
    }
    close $fh_in;
    my $guess_width = sprintf"%0d", $width_all/$count;
    return $guess_width;
}

sub usage_simple {
    die("
Usage: update_table.pl [-add|-line|-replace] [options] <STDIN|in.file>

Options: [-hit] [-m -a -b] [-d -n -A -B] [-x -i -j] [-report -change -o] 

Examples:
1. Add tabs (from -d) to file
perl update_table.pl -add -d db -A 4 -B 6 query.txt

2. Replace some tabs of file, by tabs (from -d)
perl update_table.pl -replace -a 3 -b 5 -d db -A 3 -B 5 query.txt

3. Extract records from db (-d) by ids of input
perl update_table.pl -line -m 1 -d db -n 4 query.txt 

(see -help for more details)
\n");
}

sub usage_full {
    die("
Usage: update_table.pl [-add|-line|-replace] [options] <STDIN|in.file>

Description: manipulate the tab-separated files

Options:
Type:
         -add (-D)       : if specified, add to the last column
                           ignore (-a, -b, -A, -B)
         -replace (-R)   : replace tabs by db (-d)
         -line (-L)      : extract records from db (-d) by ids of input

Common options:
         -h              : show short description
         -help           : show full description
         -hit (-U) <STR> : type of report: 0=all, 1=only hit, 2=not hit, [0]
         -d <STR>        : the db file to update input|STDIN file
         -m <INT>        : the column of [id] in db (-d), not support duplicates, [1]
         -a <INT>        : the start column in db (-d). [1]
         -b <INT>        : the end column in db (-d). [1]
         -n <INT>        : the column of [id] in input (STDIN), support duplicates, [1]
         -A <INT>        : the start column in input (STDIN), [1]
         -B <INT>        : the end column in input (STDIN), [1]

Report:
         -out (-o) <STR>    : redirect result to file
         -report (-r) <STR> : save the update type to file
         -log (-c) <STR>    : save change log to file

Extra: (optional)
         -dict <STR>    : the file to converter ids of (db) to ids of input (STDIN)
         -i <INT>       : the column of [id] db (-d) [1]
         -j <INT>       : the column of [id] input (STDIN) [2]

Examples:
1. Add tabs (from -d) to the tail of input (STDIN) file
perl update_table.pl -add -d db -A 4 -B 6 query.txt

2. Add tabs (from -d) after the col-4 of input (STDIN) file
perl update_table.pl -add -a 4 -d db -A 7 -B 12 query.txt

3. Replace some tabs of input (STDIN), by tabs (from -d) of db
perl update_table.pl -replace -a 3 -b 5 -d db -A 3 -B 5 query.txt

4. Extract records from db (-d) by ids of input (STDIN)
perl update_table.pl -line -m 1 -d db -n 4 query.txt 

5. [SORT file] Add genome positions (col-7 to 12) to query [SORT files]
perl update_table.pl -add -d db -A 7 -B 12 query.txt

6. [BED file] extract lines from db (-d) by input (STDIN, BED)
perl update_table.pl -line -m 4 -d db -n 1 query.bed

7. [use pipe] replace multiple tabs (discontinous) of input (STDIN)
perl update_table.pl -replace -m 7 -d db -A 7 -B 7 query.txt | perl update_table.pl -replace -m 9 -db -A 9 -B 9

8. use dict (-dict) option:
perl update_table.pl -add -m 4 -d db -A 7 -B 9 -dict index.txt -i 4 -j 10 query.txt

Other options:
-hit : default, 0=report all of input (STDIN)
       1=only report matched ids
       2=only report not matched ids (equal to grep -v)

Wang Ming 
wangmcas(AT)gmail.com
2015-08-08
\n");
}

__END__

Change log:

2015-07-07
  v0.1  Start this program. designed for (add, replace) purpose
        Also for the upcoming MTB db program (organize all the avaliable datasets;
        1. (Need to do) add a switch: only output the updated lines (equal extract lines from db)
        2. (Need to do) check the replaced (fields) have the same data structure.
        3. delete the criteria: -AB, -ab have to be in the same width.

2015-08-08
  v0.2
       1. error:  '-'. width to '-' x width
       2. decide the 3 purposes this script: add, replace, extract
       3. add: support add tabs to anywhere in the intput (STDIN)
       4. change: only support extract lines from db (-line), in case the format of db and input (STDIN) are differ.




### END OF FILE ###
