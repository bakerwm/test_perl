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
use Getopt::Long qw(:config no_ignore_case);
use Data::Dumper;

update_table();
exit(1);

sub update_table {
    my $info  = '';
    my $cm    = 1;
    my $cn    = 1;
    my $index = '';
    my $ci    = 1;
    my $cj    = 2;
    my $output;
    my $report_hit = 0;
#    my $not_hit;
    my $report;
    my $change;
    my $add;
    my $rp_line;
    my $replace_tabs;
    my $ca    = 1;
    my $cb    = 1;
    my $cA    = 1;
    my $cB    = 1;
    my $help;
    GetOptions(
            'f=s'          => \$info,
            'm=i'          => \$cm,
            'n=i'          => \$cn,
            'x=s'          => \$index,
            'i=i'          => \$ci,
            'j=i'          => \$cj,
            'o=s'          => \$output,
            'hit|U=i'      => \$report_hit,
#            'nothit|v'     => \$not_hit,
            'add|D'        => \$add,
            'line|L'       => \$rp_line,
            'replace|R'    => \$replace_tabs,
            'report|r=s'   => \$report,
            'change|c=s'   => \$change,
            'a=i'          => \$ca,
            'b=i'          => \$cb,
            'A=i'          => \$cA,
            'B=i'          => \$cB,
            'help|h'            => \$help,        
            ) or die "Usage: $0 [-f FILE -A 2 in.file]\n";
    usage_full() if($help);
    usage_simple() if(@ARGV == 0 && -t STDIN);

    my @inputs = ();
    while(<>) {
        chomp;
        push @inputs, $_;
    }
    # save information table to hash
    die("[-f $info] not exists\n") if(! -e $info);
    my %df = %{read_file($info, $cn)}; # save all line to values
    # convert {id} of info to {id} in <in.file> as keys.
    my %id_trans = ();
    if($index) { # information to in.file
        die("[-i,-j] Need specify -i and -j for $index\n") if(! $ci || ! $cj);
        %id_trans = %{read_file($index, $ci, $cj)}; # save one tab to values
        # update the keys of %df (to in.file ids)
        my %dn = ();
        for my $d (keys %df) {
           if(! exists $id_trans{$d}) {
next;
#               die("[$d] from (-f) not found in (-x) \n");
           }
           $dn{$id_trans{$d}} = $df{$d};
        }
        # replace the id=key in %df
        %df = ();
        %df = %dn;
    }
    # switch: which part will report: hits, not_hit, all
    die("[-hit] unknown input, expected: 0, 1 or 2\n") if($report_hit < 0 || $report_hit > 2);
    my $hit_switch = $report_hit; # 0=all, 1=only_hit, 2=not_hit
    # start conver files
    my $change_type;
    my @new = ();
    if($add) { # add info to last column
        @new = add_tabs(\@inputs, $cm, \%df, $cA, $cB, $hit_switch);
        $change_type = "Type:\tAdd columns [$cA\-$cB] to last column $info";
    }elsif($rp_line) { # replace the entire line
        @new = rp_line(\@inputs, $cm, \%df, $hit_switch);
        $change_type = "Type:\tUpdate entire lines";
    }elsif($replace_tabs) { # replace spcific columns
        @new = rp_tabs(\@inputs, $cm, \%df, $ca, $cb, $cA, $cB, $hit_switch);
        $change_type = "Type:\tReplace column [$ca\-$cb] by columns [$cA\-$cB] of ($info)";
    }else{
        die("[-add|-line|-replace] need to be specified, see -h\n");
        exit(1);
    }

# output results
    my $fh_out = *STDOUT;
    if(defined $output) {
        open $fh_out, "> $output" or die "Cannot write to $output, $!\n";
    }
    print $fh_out join("\n", @{$new[0]}) . "\n";
    close $fh_out;

# save change log to file
    my %ch_log = %{$new[1]};
    if($change) {
        open my $fh_ch, "> $change" or die "Cannot write to $change, $!\n";
        for my $c (sort {$a<=>$b} keys %ch_log) {
            print $fh_ch join("\n", $ch_log{$c}) . "\n";        
        }
        close $fh_ch;
    }

# generate a report 
    if($report) {
        open my $fh_rpt, "> $report" or die "Cannot wirte to $report, $!\n";
        my $total_lines = @inputs;
        my $update_lines = keys %ch_log;
        my $rpt = "Update:\t$update_lines of $total_lines";
        print $fh_rpt $change_type . "\n" . $rpt . "\n";
        close $fh_rpt;
    }
}

sub add_tabs {
    my @in = @{$_[0]};
    my $tm = $_[1];
    my %dd = %{$_[2]};
    my $tA = $_[3];
    my $tB = $_[4];
    my $hit_switch = $_[5];
    # 
    my $add_width = $tB - $tA + 1;
    my @add_gaps  = split //, ('-'.$add_width);
    # start a loop
    my @new_in = ();
    my %ch_log = ();
    my $count  = 1;
    for my $i (@in) {
        if($i =~ /^\#|^\s*$/) {
            push @new_in, $i;
            next;
        }
        my @ids  = extract_tabs($i, $tm, $tm);
        my $i_id = shift(@ids); 
        if(exists $dd{$i_id}) {
            next if($hit_switch == 2);
            my @add_in = extract_tabs($dd{$i_id}, $tA, $tB);
            push @new_in, join("\t", $i, @add_in);
# record the change log            
            $ch_log{$count} = "#$i\n" . join("\t", "Line-$count: Add-to-last-column:", @add_in);
        }else{
            next if($hit_switch == 1);
            push @new_in, join("\t", $i, @add_gaps);
        }
        $count ++;
    }
    return(\@new_in, \%ch_log);
}

sub rp_line {
    my @in = @{$_[0]};
    my $tm = $_[1];
    my %dc = %{$_[2]};
    my $hit_switch = $_[3];
    #
    my @tmp_val  = values %dc;
    my @dc_width = split /\t/, shift(@tmp_val);
    my @new_in   = ();
    my %ch_log   = ();
    my $count    = 1;
    for my $i (@in) {
        if($i =~ /^\#|^\s*$/) {
            push @new_in, $i;
            $count ++;
            next;
        }
        my @tabs = split /\t/, $i;
        # replace by the same width
#        die("[-f, <IN>] -f is expected not shorter than <IN>\n") if(@tabs < @dc_width);
#        die("[-f, <IN>] have different fields (width)\n") if(@tabs != @dc_width);

        my @ids  = extract_tabs($i, $tm, $tm);
        my $i_id = shift(@ids);
        if(exists $dc{$i_id}) {
            next if($hit_switch == 2);
            if(@tabs > @dc_width) {
                my $num  = @tabs - @dc_width;
                my $flag = '-' x $num;
                my @gaps = split //, $flag;
                push @new_in, join("\t", $dc{$i_id}, @gaps); # add '-' to make sure the same length
            }else {
                push @new_in, $dc{$i_id}; ## WARN: may differ in format ##
                $ch_log{$count} = "#i\nLine-$count: Update-line:\t$dc{$i_id}";
            }
        }else {
            next if($hit_switch == 1);
            push @new_in, $i;
        }
        $count ++
    }
    return(\@new_in, \%ch_log);
}

sub rp_tabs {
    my @in = @{$_[0]};
    my $tm = $_[1];
    my %df = %{$_[2]};
    my $ta = $_[3];
    my $tb = $_[4];
    my $tA = $_[5];
    my $tB = $_[6];
    my $hit_switch = $_[7];
    # check in->out in the same length
    my $in_length = $tb - $ta + 1;
    my $rp_length = $tB - $tA + 1;
#    die("[-a -b, -A -B] ($ta $tb, $tA $tB) differ in length\n") if($in_length != $rp_length);
    my @new_in = ();
    my %ch_log = ();
    my $count  = 1;
    for my $i (@in) {
        if($i =~ /^\#|^\s*$/) {
            push @new_in, $i;
            $count ++;
            next;
        }
        my @tabs = split /\t/, $i;
        my @ids  = extract_tabs($i, $tm, $tm);
        my $i_id = shift(@ids);
        if(exists $df{$i_id}) {
            next if($hit_switch == 2);
            my $pre_s = 1;
            my $pre_e = $ta - 1;
            my $nxt_s = $tb + 1;
            my $nxt_e = @tabs;
            my @pre = extract_tabs($i, $pre_s, $pre_e);
            my @nxt = extract_tabs($i, $nxt_s, $nxt_e);
            my @rep = extract_tabs($df{$i_id}, $tA, $tB);
            my @original = extract_tabs($i, $ta, $tb);
            push @new_in, join("\t", @pre, @rep, @nxt);
            $ch_log{$count} = join("\t", "#$i\nLine-$count: Replace ($ta:$tb)", @original, "by ($tA:$tB)", @rep);
        }else {
            next if($hit_switch == 1);
            next if($hit_switch);
            push @new_in, $i;
        }
        $count ++;
    }
    return(\@new_in, \%ch_log);
}

sub guess_filefmt {
# guess the input file format: PTT,SORT,BED
    my $in = $_[0];
    my $fmt = '0';
    open my $fh_in, "< $in" or die "Cannot open $in, $!\n";
    while(<$fh_in>) {
        chomp;
        next if(/^\#|^\s*$/);
        my @tabs = split /\t+/;
        next if($. < 3); # skip the first 3 lines
        last if(@tabs < 2);
        #
        if(/^\d+\.\.\d+\t[+-]\t\d+\t\d+/) {
            $fmt = 'PTT';
        }elsif(/\d+\t\d+\t[+-]/) {
            $fmt = 'SORT';            
        }elsif(/\d+\t\d+\t\w+\t\d+\t[+-]/) {
            $fmt = 'BED';
        }else {
            $fmt = '0';
        }
    }
    close $fh_in;
    return $fmt;
}

sub read_file { # file column
# save the file $in to hash, with col-$n as keys
    my ($in, $n, $k) = @_;
    $n = 1 if(! defined $n); # keys
    my %f = ();
    open my $fh_in, "< $in" or die "Cannot open $in, $!\n";
    while(<$fh_in>) {
        chomp;
        next if(/^\#|^\s*$/);
        my @tabs = split /\t/;
        die("[wrong file width or column number]\n") if($n > @tabs);
        $f{$tabs[$n - 1]} = (defined $k)?$tabs[$k - 1]:$_;
    }
    close $fh_in;
    return(\%f);
}

sub extract_tabs {
# extract the specific tabs
    my ($line, $s, $e) = @_;
    my @subs = ();
    my @tabs = split /\t/, $line;
    if($e > 0) {    
        $e = $s if($e < $s); # not pass a number to e
#        ($s, $e) = ($e, $s) if($e < $s);
        $s --;
        $e --;
        die("[wrong width of file]\n") if(@tabs < $e);
        @subs = @tabs[$s..$e];
    }
    return @subs;
}

sub replace_tabs {
# replace the a->b columns by INFO
    my ($line, $s, $e, $info) = @_;
    my @tabs = split /\t/, $line;
    $e = $s if($e < $s);
#    ($s, $e) = ($e, $s) if($s > $e);
    $s --;
#    $e;
    my @pre  = @tabs[0..$s];
    my @nxt  = @tabs[$e..$#tabs];
    return join("\t", @pre, $info, @nxt);
}

sub file_width {
# return the columns of the file
    my $in = $_[0];
    my $wide = 0;
    open my $fh_in, "< $in" or die "Cannot open $in, $!\n";
    while(<$fh_in>) {
        chomp;
        next if(/^\#|^\s*$/);
        my @tabs = split /\t/;
        if($wide != 0) {
            if($wide != @tabs) {
                die("[line $.] does not find $wide elements\n");
            }else {
                $wide = @tabs;
            }
        }
    }
    close $fh_in;
    return $wide;
}

sub usage_simple {
    die("
Usage: update_table.pl [-add|-line|-replace] [options] <STDIN|in.file>

Options: [-hit|U] [-f -n -a -b] [-m -A -B] [-x -i -j] [-report -change -o] 

1. Add tabs to the last column (col4-6 to last column)
perl update_table.pl -add -f fileB -n 1 -a 4 -b 6 fileA
perl update_table.pl -add -f fileB -n 4 -a 1 -b 6 fileA --report rpt.txt --change chg.log -o out.txt

2. update entire line
perl update_table.pl -line -f fileB -n 1 fileA 

3. Replace selected columns (update col4-6 by col4-6 of fileB)
perl update_table.pl -replace -f fileB -n 1 -A 4 -B 6 -a 4 -b 6 fileA -o out.file

(see -h for more details)
\n");
}

sub usage_full {
    die("
Usage: update_table.pl [-add|-line|-replace] [options] <STDIN|in.file>

Replace options
         -add (-D)      : if specified, add to the last column
                          ignore (-a, -b, -A, -B)
         -line (-L)     : if specified, replace the entire line of in.file by 
                          information file, ignore (-a, -b, -A, -B)
         -replace (-R)  : replace specific columns by information file

Options: <in.file>      : a file or STDIN (from pip)
         -h --help      : show this help
         -U -hit        : which part will output: 0: all, 1: only hit, 2: not hit, [0]
         -f <STR>       : file contain information for substitution
         -m <INT>       : column to search id in input file [1]
         -n <INT>       : column to search id in [-i] file [1]
         -o <STR>       : redirect output to a file
         -report (-r) <STR> : save the update results to file
         -change (-c) <STR> : save change log to file

Extra options:         
         -x <STR>       : (optional)file contain id-trans index
         -i <INT>       : column contain id for in.file [1]
         -j <INT>       : column contain id for information file [2]
         
Range of tabs to replace
         -a <INT>       : which column of input file will be repaced [1]
         -b <INT>       : (require -a) [a-b] columns will be replaced
         -A <INT>       : which column of information file for replacement [equal -a]
         -B <INT>       : (require -A) [A-B] columns for replacement [eqal -A]

[in.file] or [information] file
Tab separated file

ID-trans index
<ID-1> <ID-2>

Examples:
1. Add tabs to the last column (col4-6 to last column)
perl update_table.pl -add -f fileB -n 1 -A 4 -B 6 fileA
perl update_table.pl -add -f fileB -n 4 -A 1 -B 6 fileA --report rpt.txt --change chg.log -o out.txt
perl update_table.pl -add -f fileB -n 4 -A 1 -B 6 -x a2b.index -i 1 -j 2 -o new.fiel A.file

2. Update entire line
perl update_table.pl -line -f fileB -n 1 fileA 
perl update_table.pl -line -f fileB -n 2 -m 2 file A > new.file

3. Replace selected columns (update col4-6 by col4-6 of fileB)
perl update_table.pl -replace -f fileB -A 4 -B 6 -a 4 -b 6 fileA -o out.file

4. Update the IDs of fileA (add old ids to last column)
perl update_table.pl -add -f fileB -n 1 -A 1 -B 1 fileA | perl update_table.pl -replace -f fileB -n 1 -A 2 -B 2 -a 1 -b 1 > out.txt

Default:
(-n 1 -A 1 -B 1 -n 1 -a 1 -b 1 -i 1 -j 2)
perl update_table.pl -add -f fileB fileA     # Add columns [1-1] of fileB to last column of fileA
perl update_table.pl -line -f fileB fileA    # Update entire lines
perl update_table.pl -replace -f fileB fileA # Replace column [1-1] of fileA by columns [1-1] of fileB
\n")
}

__END__

Change log:

2015-07-07
  v0.1  Start this program. designed for (add, replace) purpose
        Also for the upcoming MTB db program (organize all the avaliable datasets;
        1. (Need to do) add a switch: only output the updated lines (equal extract lines from db)
        2. (Need to do) check the replaced (fields) have the same data structure.
        3. delete the criteria: -AB, -ab have to be in the same width.
