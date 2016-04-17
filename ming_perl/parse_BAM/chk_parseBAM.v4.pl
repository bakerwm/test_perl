#!/usr/bin/env perl

#######################################################
# parsing the BAM output from bowtie/bowtie2/tophat,  #
# 1. stat:  stat mapping reads for each sample        #
# 2. view:  create bedgraph views for each bams       #
# 3. tags:  find tags in each bam                     #
#######################################################

use strict;
use warnings;
use Cwd qw(abs_path cwd);
use File::Which;
use File::Path qw(make_path remove_tree);
use File::Basename qw(basename dirname);
use File::Spec::Functions qw(catfile catdir);
use POSIX qw(strftime);
use Getopt::Std;
use Data::Dumper;

#############################
#   >>  Main script  <<     #
#############################
&usage if(@ARGV < 1);

my $command = shift(@ARGV);
my %prog = (stat => \&statBAM,
            view => \&viewBAM,
            tags => \&BAM2tags);

die("Unknown command [$command] \n") if (!defined($prog{$command}));

#############################
# Check PERL and tools      #
# existence                 #
#############################
my %func = ("samtools"                => '',
            "bedtools"                => '',
            "htseq-count"             => '',
            "featureCounts"           => '',
            "sort2bed.pl"             => '',
            "tmp.sort2bed.pl"         => '',
            "chk_find_tags.pl"        => '',
            "search_cov_regions.pl"   => '', ## del
            "sort_to_position_sig.pl" => '',
            "update_table.pl"         => '',
            "sort2candi.pl"           => '',
            "chk_seq2rnaz.pl"         => ''
            );

my $fv = check_tools(\%func);
%func  = %{$fv};

&{$prog{$command}};

##############################
#    >>  Subroutines  <<     #
##############################

##############################
# BEGIN - 1. BAM -> mapped   #
##############################
sub statBAM {
    my %opts = (o => 'bam_stat');
    getopts('o:', \%opts);
    die(qq/
Usage: parseBAM.pl stat [options] <bam_list>

<bam_list> one BAM (full) name in one line

Options: -o    output dir, [bam_stat]
\n/) if (@ARGV == 0);

    ### Parameters
    my $bam_list  = $ARGV[0];
    my $out_dir   = abs_path($opts{o});
    my @bam_lists = @{ parse_bam_list($bam_list) };

    my $bam_stat_filename = 'bam_stat.out';
    my $bam_stat_file     = catfile($out_dir, $bam_stat_filename);
    make_path($out_dir) if( ! -d $out_dir);
    stat_bam_mapped(\@bam_lists, $out_dir, $bam_stat_file);
    print "Finish stat bam, save results to [" . $bam_stat_file . "]\n";
}

### calculate bam mapped reads
sub stat_bam_mapped {
    my @bam_lists = @{$_[0]};
    my $stat_dir  = $_[1];
    my $stat_file = catfile($stat_dir, 'bam_stat.out');
    make_path($stat_dir) if ( ! -d $stat_dir );
    open my $fh_st, "> $stat_file" or die "Cannot open file: $stat_file, $!\n";
    for my $bam (sort @bam_lists) {
        my $mapped = qx($func{'samtools'} idxstats $bam | head -n 1 | awk '{print \$3}');
        chomp($mapped);
        print $fh_st basename($bam) . "\t" . $mapped . "\n";
    }
}

### read total mapped number
sub read_bam_mapped {
    my $in = $_[0];
    my %mapped = ();
    open my $fh_in, "< $in" or die "Cannot open stat file: $in, $!\n";
    while(<$fh_in>) {
        chomp;
        next if(/^\#|^\s*$/);
        my ($bam_name, $mapped) = (split /\t/, $_)[0, 1];
        $mapped{$bam_name} = $mapped;
    }
    close $fh_in;
    return(\%mapped);    
}
##############################
# END 1. BAM -> mapped       #
##############################

##############################
# BEGIN 2. BAM -> bedgraph   #
##############################
sub viewBAM {
    my %opts = (o => 'bam_view',           
                s => 1,
                t => 'auto',
                r => 0);
    getopts('s:f:o:r:', \%opts);
    die(qq/
Usage: parseBAM.pl view [options] <bam_list>

<bam_list> one BAM (full) name in one line

Options: -o <STR>       output dir, [bam_view]
         -s <FLOAT>     scale, [0 to 1] [1]
                        0: normalize the total mapped reads to 1 million
                        1: do not normalize the mapped reads
         -f <STR>       the reference in fasta format [fa]
         -t <STR>       the bam_stat dir, [auto]
                        [auto]=redirect to the same directory of -o: [bam_stat]
                        STR=PATH\/bam_stat.out
         -r <INT>       force to calculate fwd,rev bam 0=no, 1=yes, [0]
Example:

parseBAM.pl view -o outdir -s 1 -f ref.fa inbam.list > log
\n/) if (@ARGV == 0);
    die("[-f] reference fasta file is required\n") if( ! defined $opts{f}); 
    die("[-f $opts{f}] reference file not found\n") if(! -e $opts{f});
    die("[-s $opts{s}] not support\n") if($opts{s} < 0);
    
    ### Parameters
    my $bam_list = $ARGV[0];
    my $scale    = $opts{s};
    my $ref      = $opts{f};
    my $bamstat  = $opts{t};
    my $calbam   = $opts{r};

    my $ref_idx  = abs_path($ref) . '.fai';
    my $log_filename      = '2.bam_view.log';
    my $bam_subdirname    = 'bam_subdir';
    my $bam_stat_dirname  = 'bam_stat';
    my $bam_stat_filename = 'bam_stat.out';
    my $bam_bedgraph_libname = 'bam_bedgraph.lib';

    my @bam_lists  = @{ parse_bam_list($bam_list) };
    my $out_dir    = abs_path($opts{o});
    my $log_dir    = catdir(dirname($out_dir), 'logs');
    my $log_out    = catfile($log_dir, $log_filename);
    my $bam_subdir = catdir(dirname($out_dir), $bam_subdirname);
    make_path($out_dir) if( ! -d $out_dir);
    make_path($log_dir) if( ! -d $log_dir);
    make_path($bam_subdir) if( ! -d $bam_subdir);

    ### prepare faidx file
    my $run_idx = "$func{'samtools'} faidx $opts{f}";
    if( system("$run_idx") == 0){
        ###
    }else {
        die("Failed to create ref index. \n");
    }

    ### Checking bam_stat folder
    my $bam_stat_dir  = catdir(dirname($out_dir), $bam_stat_dirname);
    my $bam_stat_file = catfile($bam_stat_dir, $bam_stat_filename);
    if($bamstat eq 'auto') {
        if( ! -e $bam_stat_file ) {
            stat_bam_mapped(\@bam_lists, $bam_stat_dir, $bam_stat_file);
        }
    }else{
        if( -e $bamstat ) {
            $bam_stat_file = $bamstat;
        }else {
            die("[-t $bamstat] file not exists, try: -t auto \n");
        }    
    }
    my %bam_mapped = %{ read_bam_mapped($bam_stat_file) };
    
    ### check bam subfolder
    ### save the fwd/rev.bam files to the directory
    my @runs = @{ check_bam_subdir($bam_subdir, \@bam_lists, $opts{r}) };
    
    ### convert bam to bedgraph files
    my $bam_bedgraph_lib = catfile($out_dir, $bam_bedgraph_libname);
    my @lib_lines = ();
    push @lib_lines, 'Generate bedgraph files using bam files';
    push @lib_lines, 'Input BAM files:';
    for my $bam (sort @bam_lists) {
        push @lib_lines, "\t" . $bam;
        push @runs, bam_to_bedgraph($bam, $ref_idx, $scale, $out_dir, \%bam_mapped);
    }
    push @lib_lines, 'Scale: 0=normalize to 1 million total reads';
    push @lib_lines, "\t" . $scale;
    open my $fh_lib, "> $bam_bedgraph_lib" or die "Cannot write to lib: $bam_bedgraph_lib, $!\n";
    print $fh_lib join("\n", @lib_lines) . "\n";
    close $fh_lib;
    ### execute commands (@runs)
    ### save the commands to $log_out file
    exe_cmd("1", \@runs);
}

### convert bam to bedgraph
sub bam_to_bedgraph {
    my $bam     = $_[0];
    my $ref_idx = $_[1];
    my $scale   = $_[2];
    my $outdir  = $_[3];
    my %bam_mapped = %{ $_[4] };
    ### get mapped reads
    my $bam_file   = basename($bam);
    my $bam_mapped = 1;
    if( exists $bam_mapped{$bam_file} ) {
        $bam_mapped = $bam_mapped{$bam_file};
    }else {
        die("[$bam] file not exists, try -t auto, and rerun this command\n");
    }
    my $sn = sprintf"%.4f", 1000000/$bam_mapped;
    $scale = $sn if($scale == 0);
    ### 
    my ($bam_name) = basename($bam) =~ /(^\w+)\./;
    my $bam_subdir = catdir(dirname($outdir), 'bam_subdir');
    my $bam_fwd = catfile($bam_subdir, $bam_name . '.fwd.bam');
    my $bam_rev = catfile($bam_subdir, $bam_name . '.rev.bam');
    my $bg_fwd  = catfile($outdir, $bam_name . '.fwd.bedgraph');
    my $bg_rev  = catfile($outdir, $bam_name . '.rev.bedgraph');
    my @runs = ();
    push @runs, "$func{bedtools} genomecov -bg -split -scale $scale -ibam $bam_fwd -g $ref_idx > $bg_fwd";
    push @runs, "$func{bedtools} genomecov -bg -split -scale $scale -ibam $bam_rev -g $ref_idx > $bg_rev";
    return(@runs);
}
##############################
# END - 2. BAM -> bedgraph   #
##############################

##############################
# BEGIN - 3. BAM -> tags     #
##############################
sub usage_bam2tags_long {
    die(qq/
Usage: parseBAM.pl tags [options] <bam_list>

<bam_list> each line contain one BAM file

Options: -o <STR>       output dir, [bam_tags]
         -f <STR>       reference file [fa]
         -g <STR>       annotation file, mtrRNA [gff]
         -n <STR>       the name of the sample [bac]
    
Type of the libraries:
         -b <STR>       libtype for each bam files. [auto]
                        auto="all libs contain >40 nt sRNA"
                        FILE="<BAM_name> <Lib_range> <Filt_range>", 
                        eg: [extent 40 nt at both ends, have to > 40 nt]
                        H37Rv01.f.s.bam 18-40  40-80
                        H37Rv02.f.s.bam 40-80  40-120
                        H37Rv03.f.s.bam 80-140 40-180
                        H37Rv04.f.s.bam >140   100-10000
         -s <FLOAT>     scale, [0 to 1] [1]
                        0: normalize the total mapped reads to 1 million
                        1: do not normalize the mapped reads
         -r <INT>       force to calculate fwd\/rev.bam? 0=no, 1=yes, [0]
         -l <INT>       force to calculate fwd\/rev.coverage? 0=no, 1=yes, [0]
         -t <STR>       bam_stat.out or "auto", [auto]
         -c <INT>       cut-off for determine edges of tags [100]
    
Control output:
         -m <INT>       merge tags from all samples, 0=no, 1=yes, [1]
         -e <INT>       count reads on each tag by htseq-count. very slow. 0=no, 1=yes, [1]
         -z <INT>       apply RNAz analysis, 0=no, 1=yes, [1]
                        Only report one region with the highest z-score ifM_name> <Lib_range> <Filt_range>,
                              272                         eg: [extent 40 nt at both ends, have to > 40 nt]
                        the input sequence contain multiple regions with 
                        z-score > 0.5. (criteria)
         -d <STR>       The database for RNAz analysis. [SixRv]
                        see: ~\/work\/database\/H37Rv\/SixRv.fa

Example:
parseBAM.pl tags -o out_dir -f ref.fa -g ref.gff inbam.list > log
\n/); 
}

sub usage_bam2tags_short {
    die(qq/
Usage: parseBAM.pl tags -o out_dir -f ref.fa -g ref.gff inbam.list > log
\n/);
}

sub BAM2tags {
    use Getopt::Long;
    my $ref;
    my $gff;
    my $out_dir       = 'bam_tags';
    my $sample_name   = 'bac';
    my $libtype       = 'auto';
    my $scale         = 1;
    my $choose_cutoff = 100;
    my $bam_stat      = 'auto';
    my $calbam        = 0;
    my $calcov        = 0;
    my $merge         = 1;
    my $exp           = 1;
    my $rnaz          = 1;
    my $rnaz_db       = '/home/wangming/work/database/H37Rv/SixRv.fa';
    my $help;
    GetOptions(
        'f|ref=s'     => \$ref,
        'g|gff=s'     => \$gff,
        'o|output=s'  => \$out_dir,
        'n|name=s'    => \$sample_name,
        'b|libtype=s' => \$libtype,
        's|scale=i'   => \$scale,
        'c|cutoff=i'  => \$choose_cutoff,
        't|bamstat=s' => \$bam_stat, # the file of bam_stat.out, or 'auot'
        'r|calbam=i'  => \$calbam,   # force to calculate fwd/rev.bam
        'l|calcov=i'  => \$calcov,   # force to calculate fwd/rev.coverage
        'm|merge=i'   => \$merge,
        'e|exp=i'     => \$exp,
        'z|rnaz=i'    => \$rnaz,
        'd|rnazdb=s'  => \$rnaz_db,
        'h|help'      => \$help,
        ) or usage_bam2tags_short();
    usage_bam2tags_long() if($help);
    usage_bam2tags_short() if(@ARGV == 0 && -t STDIN);
    die("[-f] reference fasta is required, see -h for details\n") if( ! defined $ref );
    die("[-g] gff file is required, see -g for details\n") if( ! defined $gff );
    die("[-f $ref] reference file not exist\n") if(! -e $ref);
    die("[-g $gff] annotation file not exist\n") if(! -e $gff);

    ### Parameters
    $out_dir                  = abs_path($out_dir);
    $ref                      = abs_path($ref);
    $gff                      = abs_path($gff);
    my $bam_list              = $ARGV[0];
    my @bam_lists             = @{ parse_bam_list($bam_list) };
    my $ref_idx               = $ref . ".fai";
    my $chr_name              = fetch_chrname($ref); # fetch the chrname from fasta file
    my $tag_log_out           = "3.bam_tags.log";
    my $bam_subdirname        = "bam_subdir";
    my $bam_stat_dirname      = "bam_stat";
    my $bam_stat_filename     = "bam_stat.out";
    my $bam_coverage_dirname  = "bam_coverage";
    my $bam_coverage_libname  = "bam_coverage.lib";
    my $tag_scan_dirname      = "scan_cutoff";
   
    ### Prepare working directory
    my $log_dir    = catdir(dirname($out_dir), "logs");
    my $log_out    = catfile($log_dir, $tag_log_out);
    my $bam_subdir = catdir(dirname($out_dir), $bam_subdirname);
    my $bam_coverage_dir = catdir(dirname($out_dir), $bam_coverage_dirname);
    make_path($log_dir) if( ! -d $log_dir );
    make_path($out_dir) if( ! -d $out_dir );
    make_path($bam_subdir) if( ! -d $bam_subdir );
    make_path($bam_coverage_dir) if( ! -d $bam_coverage_dir );

    ### Prepare ref.fasta index 
    my $cmd_faidx = "$func{'samtools'} faidx $ref";
    if(system("$cmd_faidx") == 0) {
        # blank !!!
    }else {
        die("Failed to create index for input ref, check: $ref\n");
    }

    ### Parse the libtype for each bam file
    my %filt_cr = %{ parse_libtype(\@bam_lists, $libtype, $out_dir) };

    ### Checking bam_stat file,
    ### generate %bam_mapped
    my $bam_stat_dir  = catdir(dirname($out_dir), $bam_stat_dirname);
    my $bam_stat_file = catfile($bam_stat_dir, $bam_stat_filename);
    if($bam_stat eq 'auto') {
        if( ! -e $bam_stat ) {
            stat_bam_mapped(\@bam_lists, $bam_stat_dir, $bam_stat_file);
        }
    }else {
        if( -e $bam_stat ) {
            $bam_stat_file = $bam_stat;
        }else {
            die("[-t $bam_stat] file not exists, or try: -t auto \n");
        }
    }
    my %bam_mapped = %{ read_bam_mapped($bam_stat_file) };

    ### Checking the bam subfolder,
    ### save the fwd/rev.bam files to the directory
    my @runs = @{ check_bam_subdir($bam_subdir, \@bam_lists, $calbam) };

    ### Checking the coverage folder,
    ### save the fwd/rev.coverage files to the directory
    my @runs_cov = @{ check_bam_coverage($bam_coverage_dir, \@bam_lists, $ref_idx, $scale, $bam_coverage_dir, $bam_coverage_libname, \%bam_mapped, $calcov) };
    push @runs, @runs_cov;
    exe_cmd("1", \@runs);

    ###############################
    # BEGIN - scan the cutoff     #
    ###############################
#    my @cutoff_lists = ('20', '50', '100', '1000', '10000', '100000');
    my @cutoff_lists = ('50', '100', '1000', '10000', '100000');
#    my @cutoff_lists = ('50', '100');
    for my $cutoff (sort {$a<=>$b} @cutoff_lists) {
        @runs = (); # clear up @runs
        ### convert bam to tags
        my $tag_scan_dir   = catdir($out_dir, $tag_scan_dirname);
        my $tag_cutoff_dir = catdir($tag_scan_dir, $cutoff);
        make_path($tag_scan_dir) if( ! -d $tag_scan_dir);
        make_path($tag_cutoff_dir) if( ! -d $tag_cutoff_dir );
        my ($run_tags, $tag_files)  = convert_bam2tag(\@bam_lists, $bam_subdir, $bam_coverage_dir, $tag_cutoff_dir, $cutoff);
        push @runs, @{$run_tags};
        exe_cmd("1", \@runs);
        ### merge/combine tags from different bam files
        # $tag_file, (ARRAY)
        combine_tags(\@{$tag_files}, $tag_cutoff_dir, $ref, $gff, \%filt_cr, $sample_name);
    }
    ###############################
    # END - scan cutoff           #
    ###############################
    
    ### Wrap files from all cutoff_folders
    ### basd on cutoff = 100
    ### $out_dir/scan_cutoff/50-100000/output/
    my $tag_type = "sRNA"; ### not used this para
    @runs = ();
    my ($run_update, $tag_update) = update_tags(\@cutoff_lists, $out_dir, $tag_type, $sample_name, $ref, $choose_cutoff);
    push @runs, @{$run_update};
    exe_cmd("1", \@runs);
    
    my @finish = ("echo \"### Finish $bam_list \"");
    exe_cmd("1", \@finish);
    ### Run exp/TPM
#    my $tag_report_dir = catdir($out_dir, 'tag_exp');
#    tag_to_count(\@bam_lists, $tag_update, );

    ### Run RNAz
#    my $tag_rnaz_dir = catdir($out_dir, 'tag_rnaz');

    
}

### count exp
#sub tag_to_count {
#    my @bam_lists  = @{ $_[0] };
#    my $tag_file   = $_[1];
#    my %bam_mapped = %{ $_[2] };
#    my $ref        = $_[3];
#
#    my $chr_name = fetch_chrname($ref);
#    for my $bam (sort @bam_lists) {
#        my ($bam_name) = basename($bam) =~ /(^\w+)/;
#        my $lib_type = "1";
#        $lib_type    = "2" if($bam_name =~ /\_[12]$/); ### featureCounts, dUTP
#        if( ! exists $bam_mapped{$bam_name} ) {
#            my $total_reads = $bam_mapped{$bam_name};
#            my $scale = sprintf"%.4f", 1000000/$total_reads;
#
#        }else {
#            die("[$bam] not found, try [chk_parseBAM.pl stat bam.list] \n");
#        }
#
#    
#    }
#
#}

### check the file libtype
### choose 'auto', or pass the criteria from a file
sub parse_libtype {
    my @bam_lists = @{ $_[0] };
    my $libtype   = $_[1]; # [auto] or a [file]
    my $out_dir   = $_[2];
    ###
    my %filt_cr = ();
    ###
    if($libtype eq "auto") {
        my $libtype_file = catfile($out_dir, "libtype.lib");
        open my $fh_out, "> $libtype_file" or die "Cannot open file $libtype_file, $!\n";
        print $fh_out join("\t", "#<BAM_name>", "<Lib_range>", "<Filt_range>") . "\n";
        for my $bam (sort @bam_lists) {
            my ($bam_name) = basename($bam) =~ /(^\w+)/;
            $filt_cr{$bam_name} = "40-100000"; # >40 nt
            print $fh_out join("\t", $bam, ">40", "40-100000") . "\n";
        }
        close $fh_out;
    }else {
        my $libtype_save = catfile($out_dir, "libtype.lib");
        open my $fh_lib, "< $libtype" or die "Cannot open [-b $libtype] file, $!\n";
        open my $fh_out, "> $libtype_save" or die "Cannot open $libtype_save, $!\n";
        while(<$fh_lib>) {
            print $fh_out $_; 
            chomp;
            next if(/^\#|\^\s*$/);
            my ($bam_file, $lib_range, $filt_range) = (split /\t/, $_);
            my ($bam_name) = basename($bam_file) =~ /(^\w+)/;
            $filt_cr{$bam_name} = $filt_range;
        }
        close $fh_lib;
    }
    ### Write a demo (H37Rv) to file, demo.libtype
    my @demo_libs = ();
    push @demo_libs, "#This is a demo file, eg: H37Rv";
    push @demo_libs, "#Header: <BAM_name>   <Lib_range> <Filt_range>";
    push @demo_libs, "#<Lib_range>: the length of RNAs for this library";
    push @demo_libs, "#<Filt_range>: filt the tags/sRNAs in this library, extent 40 nt at both ends";
    push @demo_libs, "#all tags have to larger than 40 nt";
    push @demo_libs, "#";
    push @demo_libs, "#H37Rv01.f.s.bam\t18-40\t40-80";
    push @demo_libs, "#H37Rv02.f.s.bam\t40-80\t40-120";
    push @demo_libs, "#H37Rv03.f.s.bam\t80-140\t40-180";
    push @demo_libs, "#H37Rv04.f.s.bam\t>140\t100-10000";
    push @demo_libs, "#";
    my $demo_libtype = catfile($out_dir, "libtype.demo");
    open my $fh_out, "> $demo_libtype" or die "Cannot open file $demo_libtype, $!\n";
    print $fh_out join("\n", @demo_libs) . "\n";
    close $fh_out;
    ###
    return(\%filt_cr);
}

### Check bam_coverage files
### convert the bam files to coverage files
sub check_bam_coverage {
    my $bam_covdir           = $_[0];
    my @bam_lists            = @{$_[1]};
    my $ref_idx              = $_[2];
    my $scale                = $_[3];
    my $bam_coverage_dir     = $_[4];
    my $bam_coverage_libname = $_[5];
    my %bam_mapped           = %{ $_[6] };
    my $calcov               = $_[7];
    ###
    my $bam_coverage_lib = catfile($bam_coverage_dir, $bam_coverage_libname);
    ###
    my @runs = ();
    my @lib_lines = ();
    push @lib_lines, "Generate coverage files using bam files";
    push @lib_lines, "Input BAM files:";
    for my $bam (sort @bam_lists) {
        push @lib_lines, "\t" . $bam;
        my ($bam_name) = basename($bam) =~ /(^\w+)/;
        my $bam_cov_fwd = catfile($bam_coverage_dir, $bam_name . ".fwd.coverage");
        my $bam_cov_rev = catfile($bam_coverage_dir, $bam_name . ".rev.coverage");
        ### check file existence
        next if( $calcov == 0 && -e $bam_cov_fwd && -e $bam_cov_rev);
        my @run_covs = @{ bam_to_coverage($bam, $ref_idx, $scale, $bam_coverage_dir, \%bam_mapped) };
        push @runs, @run_covs;
    }
    push @lib_lines, "Scale: 0=normaliza to 1 million total reads";
    push @lib_lines, "\t" . $scale;
    open my $fh_lib, "> $bam_coverage_lib" or die "Cannot write to lib: $bam_coverage_lib, $!\n";
    print $fh_lib join("\n", @lib_lines) . "\n";
    close $fh_lib;
    return(\@runs);
}

### convert bam to coverage
sub bam_to_coverage {
    my $bam     = $_[0];
    my $ref_idx = $_[1];
    my $scale   = $_[2];
    my $outdir  = $_[3];
    my %bam_mapped = %{ $_[4] };
    ### get mapped reads
    my $bam_file   = basename($bam);
    my $bam_mapped = 1;
    if( exists $bam_mapped{$bam_file} ) {
        $bam_mapped = $bam_mapped{$bam_file};
    }else {
        die("[$bam] file not exists, try -t auto, and rerun this command\n");
    }
    my $sn = sprintf"%.4f", 1000000/$bam_mapped;
    $scale = $sn if($scale == 0);
    ### 
    my ($bam_name) = basename($bam) =~ /(^\w+)\./;
    my $bam_subdir = catdir(dirname($outdir), "bam_subdir");
    my $bam_fwd = catfile($bam_subdir, $bam_name . ".fwd.bam");
    my $bam_rev = catfile($bam_subdir, $bam_name . ".rev.bam");
    my $bg_fwd  = catfile($outdir, $bam_name . ".fwd.coverage");
    my $bg_rev  = catfile($outdir, $bam_name . ".rev.coverage");
    my @runs = ();
    push @runs, "$func{bedtools} genomecov -d -split -scale $scale -ibam $bam_fwd -g $ref_idx > $bg_fwd";
    push @runs, "$func{bedtools} genomecov -d -split -scale $scale -ibam $bam_rev -g $ref_idx > $bg_rev";
    return(\@runs);
}

### convert bam to tags
sub convert_bam2tag {
    my @bam_lists        = @{$_[0]};
    my $bam_subdir       = $_[1];
    my $bam_coveragedir  = $_[2];
    my $tag_cutoffdir    = $_[3];
    my $cutoff           = $_[4];
    ### search tags
    my @runs = ();
    my @tag_files = ();
    for my $bam (sort @bam_lists ) {
        my ($bam_name) = basename($bam) =~ /(^\w+)\./;
        my $bam_fwd    = catfile($bam_subdir, $bam_name . ".fwd.bam");
        my $bam_rev    = catfile($bam_subdir, $bam_name . ".rev.bam");        
        my $cov_fwd    = catfile($bam_coveragedir, $bam_name . ".fwd.coverage");
        my $cov_rev    = catfile($bam_coveragedir, $bam_name . ".rev.coverage");
        my $tag_fwd    = catfile($tag_cutoffdir, $bam_name . ".fwd.tag");
        my $tag_rev    = catfile($tag_cutoffdir, $bam_name . ".rev.tag");
        my $tag_new    = catfile($tag_cutoffdir, $bam_name . ".tag.txt");
        push @runs, "perl $func{'chk_find_tags.pl'} -c $cutoff -s + $cov_fwd > $tag_fwd";
        push @runs, "perl $func{'chk_find_tags.pl'} -c $cutoff -s - $cov_rev > $tag_rev";
        push @runs, "cat $tag_fwd $tag_rev > $tag_new";
        push @runs, "# ";
        ### save tag_file to @tag_files
        push @tag_files, $tag_new;
    }
    return(\@runs, \@tag_files);
}

##############################
# ---BEGIN--- merge tags     #
##############################
### combine tags from different bam files
sub combine_tags {
    my @tag_files     = @{$_[0]};
    my $tag_cutoffdir = $_[1];
    my $ref           = $_[2];
    my $gff           = $_[3];
    my %filt_cr       = %{ $_[4] };
    my $sample_name   = $_[5];
    my @runs = ();
    my $chr_name = fetch_chrname($ref);
    ### Convert sort files to bed format.
    my ($cb_run, $tag_bed) = convert_sort2bed(\@tag_files, $chr_name);
    push @runs, @{$cb_run};

    ### Rename the tags by BAM
    my ($run_rename, $bed_renamed) = renametag($tag_bed, "1.tag_renamed");
    push @runs, @{ $run_rename };
    exe_cmd("1", \@runs);

    ### Filt the tags by criteria
    @runs = (); # clear the @runs, 
    my @tag_filts = @{ filttag($bed_renamed, "2.tag_filted", \%filt_cr) };
    push @runs, "echo \'### filt_tags in PERL script\' ";

    ### Merge multiple bed files
    my $bed_merged_dir  = catdir($tag_cutoffdir, "3.tag_merged");
    make_path($bed_merged_dir) if( ! -d $bed_merged_dir);
    my $bed_merged_file = catfile($bed_merged_dir, $sample_name . ".merged.bed");
    my @run_merge       = @{ merge_bed_files(\@tag_filts, $bed_merged_file, $ref, $gff) }; # @{$bed_renamed} !!! watch out blank bed files
    push @runs, @run_merge;
    exe_cmd("1", \@runs);

    ### Output results
    ### choose the appropriate id
    my $bed_output_dir  = catdir($tag_cutoffdir, "output");
    make_path($bed_output_dir) if( ! -d $bed_output_dir);
    report_merged_tags($bed_merged_dir, $bed_output_dir);
    
    ### report
}

### convert the txt files to bed format
sub convert_sort2bed {
    my @in_txts  = @{$_[0]};
    my $chr_name = $_[1];
    my @runs     = ();
    my @beds     = ();
    for my $in_txt (@in_txts) {
        my $out_bed = $in_txt;
        $out_bed =~ s/\.txt$/.bed/;
        push @runs, "perl $func{'tmp.sort2bed.pl'} -a sort -b bed -x 0 -n $chr_name $in_txt | sort -k 1,1 -k 2,2n > $out_bed"; # sort the bed files by "chr, start"
        push @beds, $out_bed;
    }
    return(\@runs, \@beds);
}

### add LibXX to the id of seed
sub renametag {
    my @beds         = @{ $_[0] };
    my $re_dirname   = $_[1];
    my @runs         = ();
    my @beds_renames = ();
    my $counter      = 1;
    for my $bed (@beds) {
        my $flag = sprintf"%02d", $counter;
        my $bed_rename_dir  = catdir(dirname($bed), $re_dirname);
        make_path($bed_rename_dir) if( ! -d $bed_rename_dir );
        my $bed_rename_file = catfile($bed_rename_dir, basename($bed));
        push @runs, "sed -e \'s/Seed/Lib$flag\\_Seed/\' $bed > $bed_rename_file";
        push @beds_renames, $bed_rename_file;
        $counter ++;
    }
    return(\@runs, \@beds_renames);
}

### filt the tags by criteria
### Need to execute the previous shell scripts before this step
sub filttag {
    my @in_tags      = @{ $_[0] };
    my $filt_dirname = $_[1];
    my %filt_cr      = %{ $_[2] }; # criteria 
    ### filt tags
    my @tag_filts    = ();
    for my $in_tag ( @in_tags ) {
        my $in_tag_dir    = dirname($in_tag);
        my $tag_filt_dir  = catdir(dirname($in_tag_dir), $filt_dirname);
        my $tag_filt_file = catfile($tag_filt_dir, basename($in_tag));
        make_path($tag_filt_dir) if( ! -d $tag_filt_dir);
        my ($in_tag_name) = basename($in_tag) =~ /(^\w+)/; # the name of BAM files \w+
        if(exists $filt_cr{$in_tag_name} ) {
            my ($range_low, $range_high) = (split/\:|\-/, $filt_cr{$in_tag_name});
            select_tags($in_tag, $tag_filt_file, $range_low, $range_high);
        }else {
            die("[$in_tag_name] not found in criteria text;\n");
        }
        push @tag_filts, $tag_filt_file;
    }
    return(\@tag_filts);
}

### select the tags by the criteria
### input bed format
sub select_tags {
    my $tag_in     = $_[0];
    my $tag_out    = $_[1];
    my $range_low  = $_[2];
    my $range_high = $_[3];
    my $tag_del    = $tag_out . ".del";
    ### check input: BED format
    ### check file existence
    if(not_blank_file($tag_in)) {
    unless( check_file_format($tag_in, "bed") ) {
        die("[$tag_in] file not match BED format\n");
    }
    }
    ### filt tags
    open my $fh_in, "< $tag_in" or die "Cannot open file: $tag_in, $!\n";
    open my $fh_out, "> $tag_out" or die "Cannot open file: $tag_out, $!\n";
    open my $fh_del, "> $tag_del" or die "Cannot open file: $tag_del, $!\n";
    while(<$fh_in>) {
        chomp;
        next if(/^\#|^\s*$/);
        my $length = (split /\t/, $_)[4]; # bed format: chr, start, end, id, length, strand
        if($length >= $range_low && $length <= $range_high) {
            print $fh_out $_ . "\n";
        }else {
            print $fh_del $_ . "\n";
        }    
    }
    close $fh_in;
    close $fh_out;
    close $fh_del;
}

### check the file existence, 
### BED/SORT/GFF/Fasta/blank?
sub check_file_format {
    my $file   = $_[0];
    my $format = uc($_[1]);
    ###
    ### Readin a sample of the file: the first 10 lines.
    open my $fh_in, "< $file" or die "Cannot open file: $file, $!\n";
    my $counter = 1;
    my $line    = "";
    while(<$fh_in>) {
        last if($counter > 10);
        next if(/^\#|^\s*$/);
        $line .= $_;
        $counter ++;
    }
    close $fh_in;
    ###
    ### report yes or no
    my $report = "0"; # default is no
    my $guess_fmt = "";
    ### guess the format first
    if(($line =~ s/\t/\t/g) >= 5 * ($counter - 1)) { ### tab separated file
        ### BED, SORT, GFF
        if(($line =~ /^.*\t.*\t\d+\t\d+\t\d+\t(\+|\-)/)) {
            $guess_fmt = "SORT";
        }elsif($line =~ /^.*\t\d+\t\d+\t.*\t\d+\t(\+|\-)/) {
            $guess_fmt = "BED";
        }elsif($line =~ /\d+\t\d+\t.\t(\+|\-)\t.\t/) {
            $guess_fmt = "GFF";
        }else {
            $guess_fmt = "unknown";
        }        
    }elsif($line =~ /^\>/) {
        ### fasta file
        $guess_fmt = "FASTA";
    }elsif($line eq "") {
        ### blank file
        $guess_fmt = "BLANK";

    }
    ###
#print $guess_fmt . "\n\n";
    if($guess_fmt eq $format) {
        $report = "1";
    }
    return ($report);
}

### merge the tags according to position  !!! watch out blank bed files
sub merge_bed_files {
    my @beds    = @{$_[0]};
    my $out_bed = $_[1];
    my $ref     = $_[2];
    my $gff     = $_[3];
    ### merge the files directly
    my @runs    = ();
    my $out_txt = my $out_pos = $out_bed;
    $out_txt    =~ s/\.bed/.txt/;
    $out_pos    =~ s/\.bed/.pos.txt/;
    push @runs, join(" " , 'cat', @beds, "| sort -k 1,1 -k 2,2n | $func{'bedtools'} merge -s -d -1 -c 4,5,6 -o distinct,distinct,distinct -i - > $out_bed");
#    push @runs, "perl $func{'tmp.sort2bed.pl'} -a bed -b sort $out_bed > $out_txt"; !!! not normal bed files.
    push @runs, "perl $func{'sort2bed.pl'} -t bed2sort -i $out_bed -o $out_txt";
    push @runs, "perl $func{'sort_to_position_sig.pl'} -f $ref -g $gff $out_txt > $out_pos";
    push @runs, "perl $func{'sort2candi.pl'} $out_pos";
    return(\@runs);
}

### wrap the tag files in merged folder
sub report_merged_tags {
    my $merged_dir = $_[0];
    my $report_dir = $_[1];
    ### create new id, ourput to new folder
    my @merged_files = glob("$merged_dir/*.txt");
    for my $txt (@merged_files) {
        my $txt_new = catfile($report_dir, basename($txt));
        report_tag_id($txt, $txt_new);
    }
}

### choose the first part of id
sub report_tag_id {
    my $txt_in = $_[0];
    my $txt_out = $_[1];
    open my $fh_in, "< $txt_in" or die "Cannot open file $txt_in, $!\n";
    open my $fh_out, "> $txt_out" or die "Cannot open file $txt_out, $!\n";
    while(<$fh_in>) {
        chomp;
        next if(/^\#|^\s*$/);
        my @tabs     = split /\t/, $_;
        my $id_old   = shift(@tabs);
        my ($id_new) = (split /\,|\:/, $id_old)[0];
        print $fh_out join("\t", $id_new, @tabs, $id_old) . "\n";
    }
    close $fh_in;
    close $fh_out;
}
##############################
# ---END--- merge tags       #
##############################

### update tags from each cutoff value
sub update_tags {
    my @cutoff_lists = @{ $_[0] };
    my $out_dir      = $_[1];
    my $tag_type     = $_[2]; # not used this para
    my $sample_name  = $_[3];
    my $ref          = $_[4];
    my $choose_cutoff = $_[5];
    ### Parameters
    my $chr_name        = fetch_chrname($ref);
#    my $choose_cutoff   = 100;
    my $update_dir      = catdir(catdir($out_dir, "scan_cutoff"), "update");
    my $update_cmp_dir  = catdir($update_dir, "comp");
    my $update_stay_dir = catdir($update_dir, "stay");
    my $update_out_dir  = catdir($update_dir, "output");
    make_path($update_cmp_dir) if( ! -d $update_cmp_dir);
    make_path($update_stay_dir) if( ! -d $update_stay_dir);
    make_path($update_out_dir) if( ! -d $update_out_dir);
    ### based on 100
    my $tag_filename    = $sample_name . ".merged.pos.txt";
    my $tag_choose_dir  = catdir(catdir(catdir($out_dir, "scan_cutoff"), $choose_cutoff), "output");
    my $tag_choose_file = catfile($tag_choose_dir, $tag_filename);
    my $tag_choose_bed  = $tag_choose_file;
    $tag_choose_bed     =~ s/.txt$/.bed/;
    ### runs
    my @runs = ();
    push @runs, "perl $func{'tmp.sort2bed.pl'} -a sort -b bed -x 0 $tag_choose_file > $tag_choose_bed";
    ### for each cutoff value
    my @tag_extras = ();
    for my $cutoff (sort {$a<=>$b} @cutoff_lists) {
        next if($cutoff == $choose_cutoff);
        my $tag_query_dir  = catdir(catdir(catdir($out_dir, "scan_cutoff"), $cutoff), "output");
        my $tag_query_file = catfile($tag_query_dir, $tag_filename);
        my $tag_query_bed  = $tag_query_file;
        $tag_query_bed     =~ s/.txt$/.bed/;
        ### output files
        my $cmp_bed_filename = "cutoff_" . $choose_cutoff . "_vs_" . $cutoff . ".merged.pos.bed";
        my $cmp_bed_file     = catfile($update_cmp_dir, $cmp_bed_filename);
        my $out_bed_filename = "cutoff_" . $choose_cutoff . "_vs_" . $cutoff . ".merged.pos_" . $cutoff . "_stay.txt";
        my $out_bed_file     = catfile($update_stay_dir, $out_bed_filename);
        push @tag_extras, $out_bed_file;
        ###
        push @runs, "perl $func{'tmp.sort2bed.pl'} -a sort -b bed -x 0 $tag_query_file > $tag_query_bed";
        push @runs, "$func{'bedtools'} intersect -s -wo -a $tag_choose_bed -b $tag_query_bed > $cmp_bed_file";
        push @runs, "cat $tag_query_file | perl $func{'update_table.pl'} -line -hit 2 -d $cmp_bed_file -n 10 > $out_bed_file";
    }
    ###
    my $tag_out_file = catfile($update_out_dir, $tag_filename);
    push @runs, join(" ", "cat", $tag_choose_file, @tag_extras, ">", $tag_out_file);
    push @runs, "perl $func{'sort2candi.pl'} $tag_out_file";
    return(\@runs, $tag_out_file);
}

### calculate TPM for each tags
sub txt2count {
    my $bam    = $_[0];
    my $in     = $_[1];
    my $outdir = $_[2];
    my $flag   = sprintf "%0d", $_[3];
    ### determine strand
    my $bam_name = basename($bam);
    $bam_name    =~ s/(\.|\.s.|\.f.s.|\.trim\.gz\.f\.s\.)bam//;
    my $lib_type = '1'; # SE=1, PE=2
    $lib_type    = '2' if($bam_name =~ /\_[12]$/);
    my @runs      = ();
    chomp(my $chr = qx($func{'samtools'} view $bam | head -n1 | awk '{print \$3}'));
    my $in_gff    = basename($in);
    $in_gff       =~ s/\.txt$/.gff/;
    $in_gff       = catfile($outdir, $in_gff);
    my $in_tmp    = catfile($outdir, basename($in) . '.tmp');
    my $in_tmp2   = catfile($outdir, basename($in) . '.tmp2');
    my $in_TPM    = catfile($outdir, basename($in) . '.TPM.' . $flag);
    ###
    my $fc_para   = '';
    if($bam_name  =~ /\_[1-9]$/) {
        $fc_para  = join(" ", '-M --fraction --donotsort -O -f -T 5 -g gene_id -t exon -p -P -d 40 -D 500 -s', $lib_type, '-a', $in_gff, '-o', $in_tmp);
    }else {
        $fc_para  = join(" ", '-M --fraction --donotsort -O -f -T 5 -g gene_id -t exon -s', $lib_type, '-a', $in_gff, '-o', $in_tmp);
    }
    ###
    if( not_blank_file($in) ) {
        push @runs, "perl $func{'sort2bed.pl'} -t sort2gff -f exon -s $chr -i $in -o $in_gff";
        push @runs, "$func{'featureCounts'} $fc_para $bam > $in\.log 2>&1";
        push @runs, "mv $in_tmp\.summary $in_tmp\.summary\.$flag";
        push @runs, "cat $in_tmp | sed  \'1,2 d\' | sort -k1 > $in_tmp2";
        ### count TPM
        my $bam_mapped = qx($func{'samtools'} idxstats $bam | head -n1 |awk '{print \$3}');
        my $ratio = sprintf"%.2f", 1000000/$bam_mapped;
        push @runs, "awk \'{printf(\"\%s\\t\%s\\t\%.4f\\n\", \$1, \$7, \$7*$ratio)}\'  \< $in_tmp2 | cut -f2-3 > $in_TPM";
    }
    return @runs;
}

### calculate RNAz scores for each tags
sub seq2RNAz {
    my ($txt, $ref, $rnaz_db, $outdir) = @_;
    my $txt_fa = $txt;
    $txt_fa    =~ s/\.txt$/.fa/;
    my $rnaz_log = catfile($outdir, 'rnaz.log');
#    my $rnaz_bed = catfile($outdir, 'best_RNAz.bed');
    my @runs = ();
    push @runs, "perl $func{'sort2bed.pl'} -t sort2fa -g $ref -i $txt -o $txt_fa";
    push @runs, "perl $func{'chk_seq2rnaz.pl'} RNAz -d $rnaz_db -o $outdir $txt_fa > $rnaz_log 2>&1 ";
    return @runs;
}

sub wrap_output {
    my ($info, $pos, $exp, $z) = @_;
    my $vf = fetch_id($info);
    my $vz = fetch_id($z);
    my $ve = fetch_count($exp);
    my %hf = %{$vf};
    my %hz = %{$vz};
    my %he = %{$ve};
    my $rpt_out = '';
    if( not_blank_file($pos) ) {
        open my $fh_pos, "< $pos" or die "$!";
        while(<$fh_pos>) {
            chomp;
            next if(/(^\s*$)|(^\#)/);
            my $id = (split /\t/)[0];
            my $zscore = (exists $hz{$id})?$hz{$id}:'-';
            my $note   = (exists $hf{$id})?$hf{$id}:'-';
            my $exp    = (exists $he{$id})?$he{$id}:'';
            $rpt_out  .= join("\t", $_, $exp, $zscore, $note) . "\n";
        }
        close $fh_pos;
        return $rpt_out;
    }else {
        return '';
    }
}

sub fetch_id {
    my $in   = shift(@_);
    my %info = ();
    if( not_blank_file($in) ) {
        open my $fh_in, "< $in" or die "Cannot open $in, $!\n";
        while(<$fh_in>) {
            chomp;
            my ($id, $note) = (split /\t/)[0, -1];
            if($in =~ /RNAz\.bed$/) {
                $id =~ s/^[a-zA-Z0-9]+\_//;
            }
            $info{$id} = $note;
        }
        close $fh_in;
    }
    return \%info;
}

sub fetch_count {
    my $in   = shift(@_);
    my %info = ();
    if( not_blank_file($in) ) {
        open my $fh_in, "< $in" or die "Cannot open $in, $!\n";
        while(<$fh_in>) {
            chomp;
            my @tabs = split /\t/, $_, 13;
            $info{$tabs[0]} = $tabs[-1];
        }
        close $fh_in;
    }
    return \%info;
}
##############################
# END - 3. BAM -> tags       #
##############################

##############################
# BEGIN - Global subroutines #
##############################
### Execute Shell commands
sub exe_cmd {
    my $run_switch = $_[0];
    my @runs       = @{ $_[1] };
    for my $r (@runs) {
        print $r . "\n";
        next unless($r);
        next if($r =~ /^\#|^\s*$/);
        if($run_switch) {
            if( system("$r") == 0 ){
                ##
            }else {
                warn("Failed to execute command:\n$r\n");
#                die("Failed to execute command:\n$r\n");
            }
        }
    }
}

### parse the BAM files
sub parse_bam_list {
    my $in_file = $_[0];
    die("[$in_file] bam_list file not exists\n") if( ! -e $in_file);
    my @bam_lists = ();
    open my $fh_in, "< $in_file" or die "Cannot open bam_list, $in_file, $!\n";
    while(<$fh_in>) {
        chomp;
        next if(/^\#|^\s*$/);
        die("[$_] bam file not exists in line-[$.] of $in_file\n") if( ! -e $_ );
        push @bam_lists, $_;
    }
    close $fh_in;
    return(\@bam_lists);
}

### check bam_subdir files
### split the bam files by strand, fwd/rev.bam (dUTP 2nd equal sense strand)
sub check_bam_subdir {
    my $bam_subdir = $_[0];
    my @bam_list   = @{$_[1]};
    my $cal_bam    = $_[2]; ## force to calculate fwd/rev bam
    make_path($bam_subdir) if( ! -d $bam_subdir);
    my @runs = ();
    for my $bam (sort @bam_list) {
        my ($bam_name) = basename($bam) =~ /(^\w+)\./;
        my $bam_fwd    = catfile($bam_subdir, $bam_name . '.fwd.bam');
        my $bam_rev    = catfile($bam_subdir, $bam_name . '.rev.bam');
        next if( ! $cal_bam && -e $bam_fwd  && -e $bam_rev ); ### $cal_bam = 0, ref/rev exist
        my $bam_type   = 'SE';
        $bam_type      = 'PE' if($bam_name =~ /\_[\d]$/);
        my $run        = split_bam($bam, $bam_type, $bam_fwd, $bam_rev);
        push @runs, @{$run};
    }
    return(\@runs);
}

### split BAM files by strand. (dUTP, 2nd read equal the sense strand)
sub split_bam {
    my $bam      = $_[0];
    my $bam_type = $_[1];
    my $bam_fwd  = $_[2];
    my $bam_rev  = $_[3];
    my @run_logs = ();
    if($bam_type eq 'SE') {
        push @run_logs, "$func{'samtools'} view -b -F 16 -o $bam_fwd $bam";
        push @run_logs, "$func{'samtools'} view -b -f 16 -o $bam_rev $bam";
    }elsif($bam_type eq 'PE') {
        my $temp_dir = catdir(cwd(), $$ . '_' . int(rand(1000000)));
        make_path($temp_dir);
        my $fwd1 = catfile($temp_dir, 'fwd1.bam');
        my $fwd2 = catfile($temp_dir, 'fwd2.bam');
        my $rev1 = catfile($temp_dir, 'rev1.bam');
        my $rev2 = catfile($temp_dir, 'rev2.bam');
        push @run_logs, "$func{'samtools'} view -b -f 128 -F 16 -o $fwd1 $bam";
        push @run_logs, "$func{'samtools'} view -b -f 80 -o $fwd2 $bam";
        push @run_logs, "$func{'samtools'} index $fwd1";
        push @run_logs, "$func{'samtools'} index $fwd2";
        push @run_logs, "$func{'samtools'} merge -f $bam_fwd $fwd1 $fwd2";
        push @run_logs, "$func{'samtools'} view -b -f 144 -o $rev1 $bam";
        push @run_logs, "$func{'samtools'} view -b -f 64 -F 16 -o $rev2 $bam";
        push @run_logs, "$func{'samtools'} index $rev1";
        push @run_logs, "$func{'samtools'} index $rev2";
        push @run_logs, "$func{'samtools'} merge -f $bam_rev $rev1 $rev2";
        push @run_logs, "rm -r $temp_dir";
    }else{
        ### 
    }
    return(\@run_logs);
}

### fetch the chromosome name from fasta file
sub fetch_chrname {
    my $in_fa = $_[0];
    my $chr_name = '';
    open my $fh_in, "< $in_fa" or die "Cannot open file: $in_fa, $!\n";
    while(<$fh_in>) {
        chomp;
        if(/^\>/) {
            if(/^\>gi\|/) {
                $chr_name = (split /\|/, $_)[3]; # choose NC_123456.1
            }elsif(/\s/) {
                $chr_name = (split /\s/, $_)[0]; # choose the first part of the id
            }else{
                $chr_name = $_;
            }
            last;
        }
    }
    close $fh_in;
    $chr_name =~ s/^\>//; # remove the header flag
    return($chr_name);
}

### check file blank or not
sub not_blank_file {
    my $in = shift(@_);
    if( -e $in ) {
        open my $fh_in, "< $in" or die "$!";
        my $count = 0;
        while(<$fh_in>) {
            chomp;
            next if(/(^\s*$)|(^\#)/); # skip blank or comment lines
            $count ++;
        }
        close $fh_in;
        return $count;
    }else {
        return 0;
    }
}

### Search the tools, perl scripts in the following PATH:
### 1. ENV (in the default ENV, such as, bedtools, samtools, RNAz)
### 2. perl scripts in '~/work/bin/' and '~/work/bin/temp/'
sub check_tools {
    my %func = %{ $_[0] };
    # @_ input the name of tool
    my @missing;
    for my $t (sort keys %func) {
        if( tool_path($t) ) {
            $func{$t} = tool_path($t);
        }else {
            push @missing, $t;
        }
    }
    if(@missing) {
        printf("The following tools are missing \n");
        for my $p (sort @missing) {
            printf("%-15s : Not found in \$PATH, %-30s\n\n", $p, '~/work/bin/temp/');
        }
    }
    my $st = (@missing)?0:1;
    return ($st, \%func);
#
    sub tool_path {
        my $tool     = $_[0];
        my $perldir  = $ENV{HOME} . '/work/bin/temp';
        my $perldir2 = $ENV{HOME} . '/work/bin';
        if(-e catfile($perldir, $tool)) {
            return catfile($perldir, $tool);
        }elsif( -e catfile($perldir2, $tool) ) {
            return catfile($perldir2, $tool);
        }elsif(which($tool)) {
            return $tool;
        }else {
            return 0;
        }
    }
}

sub usage {
    die(qq/
Usage: parseBAM.pl <command> [arguments]\n
Command: stat   count mapping reads for each BAM
         view   create "*.bedgraph" files for each BAM
         tags   find tags from each BAM files (sRNA candidates)
\n/);
}
##############################
# END - Global subroutines   #
##############################

# Structure of the output directory:
# Output
#   |-bam_stat
#   |   |-bam_stat.out (output - 1)
#   |
#   |-bam_subdir
#   |   |-*.fwd/rev.bam (bam files)
#   |
#   |-bam_coverage
#   |   |-*.fwd/rev.coverage
#   |
#   |-bam_view
#   |   |- *.fwd.bedgraph, *.rev.bedgraph, *.bam (output - 2)
#   |
#   |-bam_tags
#   |   |-libtype.lib (library type for each bam files)
#   |   |-scan_cutoff
#   |   |   |-50
#   |   |   |-100
#   |   |   |-1000
#   |   |   |-10000
#   |   |   |-100000
#   |   |   |-update
#   |   |   |   |-comp
#   |   |   |   |-stay
#   |   |   |   |-output (contain tag files)
#   |   
#   |-report
#   |   |-Lib01.report.txt, Lib02.report.txt, ...       
#
# change log
#
# 2012-05-01
# 1. convert bedgraph / hitogram graph using 'bedtools genomecov'
# 2. for strand-specific RNA-Seq, split PE BAM into two files: fwd.bam and rev.bam
# 3. find cov-tags with para: -cut-off: 100,
# 4. Using 6 MTB Complex genomes for RNAz analysis:
#      NC_008769: Mycobacterium bovis BCG str. Pasteur 1173P2
#      NC_015848: Mycobacterium canettii CIPT 140010059
#      NC_008596: Mycobacterium smegmatis str. MC2 155
#      NC_009525: Mycobacterium tuberculosis H37Ra
#      NC_000962: Mycobacterium tuberculosis H37Rv
#      NC_012943: Mycobacterium tuberculosis KZN 1435
# 5. set RNAz --cut-off=0.5
# 6. split these program into 3, for different purpose
#      stat: calculate total mapped read
#      view: create bedgraph files for genome browsers (eg: artemis, IGB)
#      tags: find cov regions in BAM files.
# 7. RNAz output: only report one region with the highest z-score if the input seq 
#    contain multiple regions with z-score > 0.5.
# 8. find cov regions (sRNA candidates): candidate should be 60 bp and 100 bp to its 
#    neighbor CDSs (genes).
# 9. you can input a custom GFF files (eg: only contain CDSs), to find sRNA candidates.
# 10. wrap count, TPM and RNAz score to one file in report directory.
# 
# 2014-11-03
# v0.1
#     1. support multiple bam files 
#     2. ONLY support single-chromosome sample
#     3. merge_tags by bedtools mrege: -s -d -1 -c 4,5,6 -o distinct,distinct,distinct
#     4. output results in out.dir/03.seqs
#
# 2015-02-07
# v0.2
#     1. support single-bam input file, Named: LibN
#     2. delete the step: copy reference file to current dir. insteat read the original fa/gff files
#
# 2015-04-05
# v0.3
#     1. Using 'HTSeq-count' instead of "bedtools multicov" to count reads on each features
#     2. Splite the BAM file into strand-specific files: fwd.bam and rev.bam
#        SE: 
#            samtools view -b -F 16 in.bam -o fwd.bam
#            samtools view -b -f 16 in.bam -o rev.bam
#        PE:
#            samtools view -b -f 128 -F 16 in.bam -o fwd1.bam
#            samtools view -b -f 80 in.bam -o fwd2.bam
#            samtools index fwd1.bam
#            samtools index fwd2.bam
#            samtools merge fwd.bam fwd1.bam fwd2.bam
#            #
#            samtools view -b -f 144 in.bam -o rev1.bam
#            samtools view -b -f 64 -F 16 in.bam -o rev2.bam
#            samtools index rev1.bam
#            samtools index rev2.bam
#            samtools merge rev.bam rev1.bam rev2.bam
#    3. Perform multiple sequence alignment by: ClustalW2 version 2.1, with default parameter
#    4. Perform RNAz analysis (RNAz version 2.1) --cut-off=0.5
#
# 2015-06-02
# v0.4
#    1. count reads on features by HTSeq-count: (find online note:)
#        SE: --stranded=yes
#        PE: --stranded=reverse (dUTP ssRNA-Seq)
#
# bug:
#    Line188: if input txt is blank, skip the 'paste' step, and copy the _sRNA to _count file. [BUG, not fix ye]
#
# 2015-06-15
# v0.5
#    1. replace HTSeq-count by FeatureCounts to count reads on each sRNAs. (much faster)
#    2. change para for featureCounts: -M --fraction count reads 1/n
# 
# 2015-12-30
# v0.6
#    1. Reconstruct this perl script, delete the time consuming steps: generating fwd/rev.bam, fwd/rev.coverage files.
#    2. Introduce a new file_folder structure: "bam_stat", "bam_subdir", "bam_coverage", "bam_view", "bam_tags"
#    3. Scan the cutoff values in this perl script
#    4. Support input lib_file that contain the library types for each bam file
#
# Author: Wang Ming, wangmcas(AT)gmail.com
### END ###
