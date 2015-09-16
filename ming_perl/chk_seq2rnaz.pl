#!/usr/bin/env perl
#
################################################
# fa to RNAz analysis
# 1. Extract the best hits from database 
# (keep only one hit for each input seq)
# 2. Perform multiple sequences alignment using 
# ClustalW2, default parameter
# 3. RNAz analysis using the following parameters:
#   - rnazWindow.pl : --min-seqs=2 --max-seqs=6
#   - RNAz : --forward --no-shuffle --cutoff=0.5
#   - rnazCluster.pl, rnazIndex.pl, rnazBEDsort.pl 
#     default parameters
# 4. wrap all output to "best_RNAz.bed"
#
# Wang Ming wangmcas(at)gmail.com
#################################################

use strict;
use warnings;

use Cwd;
use File::Which;
use File::Path qw(make_path remove_tree);
use File::Spec::Functions qw(catdir catfile);
use File::Basename qw(basename dirname);
use Bio::SearchIO;
use Bio::SeqIO;
use POSIX qw(strftime);
use Getopt::Std;
use Data::Dumper;

my @modules = ("Cwd",
               "File::Which", 
               "File::Path qw(make_path remove_tree)",
               "File::Spec::Functions qw(catdir catfile)",
               "File::Basename qw(basename dirname)",
               "Bio::SearchIO",
               "Bio::SeqIO",
               "POSIX qw(strftime)",
               "Getopt::Std");
 
my %func = ("sort2bed.pl" => '', 
            "blastall"    => '',
            "formatdb"    => '',
            "clustalw2"   => '',
            "RNAz"        => '',
            "rnazWindow.pl"  => '',
            "rnazCluster.pl" => '',
            "rnazIndex.pl"   => '',
            "rnazBEDsort.pl" => '');

&usage if(@ARGV < 1);

my $command = shift(@ARGV);
my %prog    = (RNAz   => \&runRNAz, 
               module => \&checkModules,
               demo   => \&runDemo);
die("Unkown command [$command]\n") if(!defined($prog{$command}));

&{$prog{$command}};
#exit(0);

###
sub checkModules {
    my $flag = shift(@_);
    my ($fst, $fha) = check_tools(\%func);
    my ($pst, $line)  = check_modules('modules');
    my $funclist = '';
    my $check_log = '';
    %func = %{$fha};
    for my $f (sort keys %func) {
        $funclist .= sprintf("%-15s : %-40s\n", $f, $func{$f});
    }
    if($fst && $pst) {
        $check_log .= "Everything is OK\n";
        $check_log .= "[Perl Modules]\n";
        $check_log .= $line. "\n";
        $check_log .= "[Perl scripts and Commands]\n";
        $check_log .= $funclist. "\n";  
    }else {
        $check_log .= "Something wrong:\n";
        $check_log .= "Perl modules status:\n\n";
        $check_log .= $line. "\n";
        $check_log .= "Command or Perl script status:\n\n";
        $check_log .= $funclist. "\n\n";
    }
    if($flag) {
        return $check_log;
    }else {
        print $check_log;
    }
}
    
###
sub runRNAz {
    my %opts = (d => '/home/wangming/work/database/H37Rv/SixRv.fa', o => 'RNAz_out');
    getopts("d:o:", \%opts);
    die(qq/
Usage: seq2RNAz.pl [options] <input.fa>

Options: -d :   The database file for multiple sequences alignment. [sixMTB]
         -o :   The output directory, [RNAz_out]

Note: For the input sequences that contain multiple regions meet the RNAz criteria, only the
      region with the largest RNAz score is reported. [out_path\/best_RNAz.bed]
\n/) if(@ARGV == 0);
    
    my $tmplog = checkModules(1);
    die("[-d] file not exists or not a fasta\n") if(! (-e $opts{d} && check_fasta($opts{d}, 0)) );
    make_path($opts{o}) if(! -d $opts{o});
#    die("[-o] require the output dir\n") if(! -d $opts{o});
# clean the SeqFA dir.
    remove_tree(catdir($opts{o}, 'SeqFA'));
    my $infa     = shift(@ARGV);
    my $hit_seq  = run_blast($infa, $opts{d}, $opts{o});
    my $rnaz_out = '';
    my $rnaz_file = catfile($opts{o}, 'best_RNAz.bed');
    for my $h (sort keys %{$hit_seq}) {
        my $hseq = catfile(catdir($opts{o}, 'SeqFA'), $h.'.fa');
        my $haln = catfile(catdir($opts{o}, 'SeqFA'), $h.'.aln');
        my $zbed = catfile(catdir($opts{o}, 'SeqFA'), $h.'.bed');
        if( run_clustalw2($hseq) ) {
           if( run_rnaz($haln) ) {
               if( my $zout = wrap_rnazout($zbed) ){
                   $rnaz_out .= $zout;
               }
           }
        }
    }
    open my $fh_z, "> $rnaz_file", or die "Cannot open $rnaz_file, $!\n";
    print $fh_z $rnaz_out;
    close $fh_z;
}

### 
sub runDemo {
    my %opts = ();
    getopts("o:", \%opts);
    die "[-o] Need input a output dir:\n" if (! defined $opts{o});
    make_path($opts{o}) if(! -d $opts{o});
    my $tmp_log = checkModules(1);
    # test seq from H37Rv (NC_000962.3)
    # 1. Rvnr03 is one of the ribosome RNAs in H37Rv. [P]
    # 2. Rvnt01 is one of the tRNAs in H37Rv. [P]
    # 3. rand_1074903_100, is a random selected 100-nt sequence in H37Rv (start at : 10749903). [N]
    my $db = '/home/wangming/work/database/H37Rv/SixRv.fa';
    my $demo_fa = '>Rvnr03
TTACGGCGGCCACAGCGGCAGGGAAACGCCCGGTCCCATTCCGAACCCGGAAGCTAAGCCTGCCAGCGCCGATGATACTGCCCCTCCGGGTGGAAAAGTAGGACACCGCCGAACA
>Rvnt01
GGGCCTATAGCTCAGGCGGTTAGAGCGCTTCGCTGATAACGAAGAGGTCGGAGGTTCGAGTCCTCCTAGGCCCA
>rand_1074903_100
CGATCATCGCCCTGATGGTGGCGTCGAGGTTGGCGAGCTGCTGCTGGACTGTCTCAAGATCGGGTCGTCCGTTGACGATTTTCTGCCGGCGGTCCAGCTC';
    my $rnaz_dir = $opts{o};
    my $infa     = catfile($rnaz_dir, 'input.fa');
    open my $fh_in, "> $infa" or die "Cannot open $infa, $!\n";
    print $fh_in $demo_fa . "\n";
    close $fh_in;
    my $hit_seq  = run_blast($infa, $db, $rnaz_dir);
    my $rnaz_out = '';
    for my $h (sort keys %{$hit_seq}) {
        my $hseq = catfile(catdir($rnaz_dir, 'SeqFA'), $h.'.fa');
        my $haln = catfile(catdir($rnaz_dir, 'SeqFA'), $h.'.aln');
        my $zbed = catfile(catdir($rnaz_dir, 'SeqFA'), $h.'.bed');
        if( run_clustalw2($hseq) ) {
            if( run_rnaz($haln) ) {
                if( my $zout = wrap_rnazout($zbed) ){
                    $rnaz_out .= $zout;
                }
            }
        }
    }
    print '[Input] 3 sequences:'. "\n";
    print $infa . "\n\n";
    print '[Expect] '. "\n". 'Only the first 2 sequences (Rvnt01 and Rnvr01) will print out:' . "\n\n";
    print '[Result]' . "\n" . $rnaz_out . "\n";
}

sub run_blast {
    my ($in, $db, $outdir) = @_;
    my $out_blast = catfile($outdir, 'out.blast');
    system "$func{'formatdb'} -p F -i $db";
    system "$func{'blastall'} -p blastn -i $in -d $db -e 1e-2 -m 7 -o $out_blast";
    my $seq_num = parse_blast($out_blast);
    return $seq_num;
}

sub parse_blast{ # extract one best hit for each genome
    my $blast_out = shift(@_); # xml format
    my $in = Bio::SearchIO->new(-format => 'blastxml',
                                -file => $blast_out);
    my $best_hits = catfile(dirname($blast_out), 'best_hits.txt'); ## output best hits txta
    my $seq_dir = catdir(dirname($blast_out), 'SeqFA');
    make_path($seq_dir) if(! -e $seq_dir);
    my %seq_count;
    open my $fh_hits, "> $best_hits" or die "$!";
    while (my $result = $in->next_result) {
        next if($result->num_hits == 0);
        my $q_name = $result->query_name;
        my $q_length = $result->query_length;
        my $q_acc = $result->query_accession;
        my $count = 0;
        while(my $hit = $result->next_hit) {
            my $s_name = $hit->name;
            my $s_length = $hit->length;
            my $s_sig = $hit->significance;
            my $s_description = $hit->description;
            my $hsp_num = $hit->num_hsps;
            # Choose the first HSP hit;
            my $hsp = $hit->next_hsp;
            my ($q_start, $q_end) = ($hsp->start('query'), $hsp->end('query'));
            my ($s_start, $s_end) = ($hsp->start('hit'), $hsp->end('hit'));
            my $s_evalue = $hsp->evalue;
            my $s_bits = $hsp->bits;
            my $frac_identical = sprintf "%.1f", $hsp->frac_identical * 100;
            my $frac_conserved = sprintf "%.1f", $hsp->frac_conserved * 100;
            my $hsp_string = $hsp->hit_string;
            if($hsp->strand('hit') eq '-1') {
                $hsp_string = revcomp($hsp->hit_string);
                ($s_start, $s_end) = ($s_end, $s_start);
            }
            if($hsp->strand('query') eq '-1') {
                ($q_start, $q_end) = ($q_end, $q_start);
            }
            # Create new strain name;
            my ($strain) = $s_description =~ /\w+\s\w+\s(\w+)/;
#            $strain = 'Msmeg' if ($strain eq 'str');
            my $s_ID = $strain.'_'.$q_name; # the strain name need to be trimed
            ## output fasta files
            $seq_count{$q_name} ++;
            my $seq_fa = catfile($seq_dir, "$q_name\.fa"); # multiple seqs
            open my $fh_fa, ">> $seq_fa" or die "$!" ;
            print $fh_fa "\>$s_ID\n$hsp_string\n";
            close $fh_fa;
            my $q_info = join "\t", ($q_name, $s_name, $q_length, $s_length, $frac_identical,
                                     $q_start, $q_end, $s_start, $s_end, $s_evalue, $s_bits);
            print $fh_hits $q_info, "\n";
        }
    }
    close $fh_hits;
#    return $best_hits_file;
    return (\%seq_count);
}

sub revcomp {
    my $in = shift(@_);
    my $out = $in;
    $out =~ tr/ATCGatcg/TAGCtagc/;
    $out = reverse($out);
    return $out;
}

sub run_clustalw2 {
    my $in_fa  = shift(@_);
    if( check_fasta($in_fa, 1) ) {
        system"mv -f $in_fa\.tmp $in_fa";
        my $clustal_log = $in_fa . '.clu.log'; #catfile(dirname($in_fa), $in_fa . 'clu.log');
        system "$func{'clustalw2'} -INFILE=$in_fa > $clustal_log 2>&1 ";
print "$func{'clustalw2'} -INFILE=$in_fa > $clustal_log 2>&1 " . "\n";
        return 1;
    }else {
        return 0;
    }
}

sub check_fasta { 
# warn, manipunate original file
# 1. delete the short sequence: < 10 nt
# 2. delete the blank sequence
# 3. return the number of seqs
# 4. wirte the seq to new file "in.fa.tmp"
    my ($seq, $type) = @_;
    my $seq_tmp = $seq. ".tmp";
    my $fa_out_num = 0;
    my $seqin = Bio::SeqIO->new(-file    => "< $seq",
                                -format  => 'fasta');
    my $seqout = Bio::SeqIO->new(-file   => "> $seq_tmp",
                                -format  => 'fasta');
    while(my $s = $seqin->next_seq) {
        my $newid = (split /\,|\:/, $s->id)[0]; # the first id
        $s->desc(''); # delete description field
        if(length($s->seq) > 10 ) {
            $fa_out_num ++;
            if($type) {
                $seqout->write_seq($s);
            }
        }
    }
    return $fa_out_num;
}

sub run_rnaz {
    my $in_aln  = shift(@_);
    my $out_bed = catfile(dirname($in_aln), basename($in_aln));
    $out_bed    =~ s/\.aln/.bed/;
    my $cmd_line = join(" ", $func{'rnazWindow.pl'}, "--min-seqs=2 --max-seqs=6",
                        $in_aln, "| RNAz --forward --no-shuffle ", # --cutoff=0.5
                        "| rnazCluster.pl | rnazIndex.pl --bed | rnazBEDsort.pl",
                        ">", $out_bed);
    system "$cmd_line";
    return  1;
}

sub wrap_rnazout {
# input the bed file of RNAz result
# 1. choose the best z-score for each input seq
# 2. return the output
    my $in = shift(@_);
    my %hit = ();
    open my $fh_in, "< $in" or die "Cannot open $in, $!\n";
    while(<$fh_in>){
        chomp;
        s/^\>//;
        next if (/^\s*$/); # skip blank lines
        my @tabs = split /\t/;
        my $id   = shift(@tabs);
        $id = (split /\_/, $id, 2)[1]; # trim the prefix of IDs, eg: KZN_ncRv0001 => ncRv0001
        push @{$hit{$id}}, join("\t", $id, @tabs);
    }
    close $fh_in;
    my $out = '';
    if((keys %hit)) {
        my ($best, $best_des) = (0.5, '');
        for my $n (keys %hit) {
            for my $f (@{$hit{$n}}) {
                my $score = (split /\t/, $f)[-1];
                if($score > $best) {
                    ($best, $best_des) = ($score, $f);
                }
            }
            $out .= $best_des . "\n";
        }
        return $out;
    }else {
        return 0;
    }
}

sub check_tools {
    my $tool = shift(@_);
    my %func = %{$tool};
    # @_ input the name of tool
    my @missing;
    for my $t (keys %func) {
        if( tool_path($t) ) {
            $func{$t} = tool_path($t);
        }else {
            push @missing, $t;
        }
    }
    if(@missing) {
        printf("The following tools are missing in both \$PATH and ~/work/bin/temp\n");
        for my $p (sort @missing) {
            printf("%-15s : Not found in \$PATH, %-30s\n\n", $p, '~/work/bin/temp/');
        }
    }
    my $st = (@missing)?0:1;
    return ($st, \%func);
#
    sub tool_path {
        my $tool = shift(@_);
        my $perldir = $ENV{HOME}. '/work/bin/temp';
        if(-e catfile($perldir, $tool)) {
            return catfile($perldir, $tool);
        }elsif(which($tool)) {
            return $tool;
        }else {
            return 0;
        }
    }
}

sub check_modules {
    my $list = '';
    my @missing;
    for my $m (@modules) {
        my ($root) = split(" ", $m);
        if(eval "use $m; 1") {
            #
        } else {
            push @missing, [$m, $@];
        }
    }
    if(grep($_ =~ /modules?/i, @_)) {
        for my $m (sort @modules) {
            my $is_missing = grep($_->[0] eq $m, @missing) || 0;
            $m =~ s/ .*//;
            my $version = "";
            if(! $is_missing) {
                eval {my $ver_str = $m . "::VERSION"; $version = eval "\$$ver_str"};
                $version ||= "?";
            }
            $list .= sprintf("%-8s %8s %s\n", $is_missing ? "missing" : "ok", $version, $m);
        }
    }
    if(@missing) {
        $list .= sprintf("\n*** Required module(s) missing ***\n\n");
        $list .= sprintf("You can install the modules by \n\$ cpan File::Which \nor \n\$ cpanm File::Which \n\nin your console\n\n");
        for my $pair (sort @missing) {
            my ($m, $error) = @$pair;
            $m =~ s/ .*//;
            $list .= sprintf("missing %s\n", $m);
            $error =~ s/\n//g;
            $error =~ s/BEGIN failed.*//;
            $list .= sprintf(" error %s\n", $error);
        }
    }
    my $st = (@missing)?0:1;
    return ($st, $list);
}

sub usage {
    die(qq/
Usage: seq2rnaz.pl <command> [arguments]

Command: RNAz       Perform RNAz analysis for the input seqeuence
         demo       Run a demo to test the program
         module     Checking the required Perl Modules in PATH
\n/);
}

__END__

Change log
1. skip the blank lines in RNAz report
2. skip the short seq (<10 nt) in input fa, and error fasta files
3. correct the reverse strand seq.
4. add the step: check required Perl modules and perl scripts
5. add a demo, 3 (2 RNAs and 1 random seq) seqs of H37Rv to test the program
6. wrap the RNAz output (bed), skip blank lines.
7. only report the region with the largest RNAz score, if multiple
   regions of the input fa match the criteria: RNAz score > 0.5.

2015-07-09
   1. add "-INFILE=in.fa" for clustalw2 program. (error: not recognize the input type)

2015-08-08
   1. trim the prefix (strain name) in the output of RNAz.
   eg: KZN_Rv0001 => Rv0001


