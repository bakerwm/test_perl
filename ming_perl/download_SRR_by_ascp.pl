#!/usr/bin/perl -w

##################################################
# Download NCBI data using Aspera                #
#                                                #
# 1. SRA data.                                   #
# 2. Bacteria genomes.                           #
#                                                #
# Wang Ming wangmcas(AT)gmail.com                #
# 2015-06-22                                     #
##################################################

use strict;
use warnings;
use File::Basename qw(basename dirname);
use File::Spec::Functions qw(catdir);
use File::Path qw(make_path);
use Getopt::Std;

my $ascp_exe = '~/.aspera/connect/bin/ascp';
my $ascp_key_openssh = '~/.aspera/connect/etc/asperaweb_id_dsa.openssh';

ascp_download();
exit(1);

# main script
sub ascp_download {
    my %opts = (t => 'bac', m => '1'); # tyep: sra, bac
    getopts("t:o:m:d:", \%opts);
    usage() if(@ARGV != 1);
    die("[-o] Need specify the output dir:\n") if(! defined $opts{o});
    die("[-t $opts{t}] Unknown type. support: [sra, bac]") if(! $opts{t} =~ /^(sra|bac)$/i);
    make_path($opts{o}) if(! -d $opts{o});
    # search mode
    die("[-m $opts{m}] Unknown mode. support: [1, 2]\n") if(! $opts{m} =~ /^(1|2)$/);
    my $input = shift(@ARGV);
    my @ids   = read_input($input);
    if($opts{m} eq '2') {
        die("[-d] Need input the index of bacteri genome info\n") if(! defined $opts{d});
        @ids = search_ids(\@ids, $opts{d});
    }
    # run ascp downlooad program:
    print 'Found '. scalar(@ids) . ' ids:' . "\n";
    my @runs  = id_to_run(\@ids, $opts{t}, $opts{o});
    for(my $i = 0; $i <= $#ids; $i ++) {
        my $count = $i + 1;
        print '['. $count . '/'. scalar(@ids) .'] ' .  $ids[$i]. "\n";
#        print $runs[$i]. "\n";
        system "$runs[$i]";
    }
}

# search genome names by input queries
# eg: H37Rv, Mycobacterium, ...
sub search_ids {
    my ($query, $index) = @_;
    my @hit_genomes = ();
    my %info = ();
    open my $fh_in, "< $index" or die "Cannot open $index, $!\n";
    while(<$fh_in>) {
        chomp;
        next if(/(^\s*$)|^\#/);
        my @tabs = split /\s+/, $_;
        my $note = '-';
        $note = $tabs[1] if (defined $tabs[1]);
        $info{$tabs[0]} = $note;
    }
    close $fh_in;
    # search id
    for my $n (@$query) {
        my @hits_1 = grep { index($_, $n) > -1 } keys(%info);
        my @hits_2 = grep { index($info{$_}, $n) > -1 } keys(%info);
        push @hit_genomes, @hits_1, @hits_2;    
    }
    return @hit_genomes;
}

# read input para: id list
sub read_input {
    my $in = shift(@_);
    my %ids = ();
    if( -e $in ){
        open my $fh_in, "< $in" or die "$!\n";
        while(<$fh_in>) {
            chomp;
            next if(/(^\#)|(^\s*$)/);
            my $id = (split /\s+/, $_)[0];
            $ids{$id} ++;
        }
        close $fh_in;
        die("[$in] : No ids found\n") if((keys %ids) < 1);
    }else {
        $ids{$in} ++;
#        push @ids, $in;
    }
    return (sort keys %ids);
}

# generate ascp commands for all input ids
sub id_to_run {
    my ($id, $t, $out_path) = @_;
    my @runs    = ();
    for my $i (sort @$id) {
        my $r = id_to_ascprun($i, $t, $out_path);
        push @runs, $r;
    }
    return @runs;
}

# generate command for one id
sub id_to_ascprun {
    my ($id, $type, $out_path) = (@_);
    my $ncbi_ascp_url = '';
    if($type eq 'sra') {
        $ncbi_ascp_url = srr_id_url($id);
    }elsif($type eq 'bac') {
        $ncbi_ascp_url = bac_genome_url($id);
#        $out_path = catdir($out_path, $id);
#        make_path($out_path) if(! -d $out_path);
    }else {
        #
    }
    my $ncbi_ascp_para   = join(' ', $ascp_exe, '-i', $ascp_key_openssh, '-q -k 1 -T -l200m');
    my $ncbi_ascp_domain = 'anonftp@ftp.ncbi.nlm.nih.gov:';
    my $cmd_line = join(' ', $ncbi_ascp_para, $ncbi_ascp_domain . $ncbi_ascp_url, $out_path);
    return $cmd_line;
}

# generate url for SRA data
sub srr_id_url {
    my $id = shift(@_);
    my $srr_path_prefix  = '/sra/sra-instant/reads/ByRun/sra';
    my ($pre_1, $pre_2)  = $id =~ /(SRR|ERR|DRR)(\d{3})/;
    my $ncbi_srr_url    = join('/', $srr_path_prefix, $pre_1, $pre_1 . $pre_2, $id, $id.'.sra');
    return $ncbi_srr_url;
}

# generate url for bacteria genomes
sub bac_genome_url { # download the folder
    my $id = shift(@_);
    my $bac_path_prefix  = '/genomes/Bacteria';
    my $ncbi_bac_url     = join('/', $bac_path_prefix, $id);
    return $ncbi_bac_url;
}

sub usage {
    die("
Usage: ascp_download.pl [options] <input>

Options: -o <str>   : output dir for download files
         -t <str>   : Type of download file: 
                      sra=SRA data, bac=genome data, [default: bac]
         -m <int>   : the mode to search the string name, 
                      1=full name, 2=part name, [default: 1]
         -d <str>   : the index file contain full bacteria list in ncbi
                      see: [ncbi_bac_genome_list.pl]
         <input>    : input can be a id, or a file contain a list of ids (one in each line);
#
Example:
1. download SRA file by the SRA id.
ascp_download.pl -t sra -o out SRR123456

2. download bacteria genomes by the name.
ascp_download.pl -t bac -o out Mycobacterium_JDM601_uid67369 

3. search bacteria genomes, and download.
ascp_download.pl -t bac -o out -d bacteria.info -m search  H37Rv ...
\n");
}

