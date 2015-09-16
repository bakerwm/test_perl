#!/usr/bin/perl -w

#################################################
# Find GO annotation of input genes
#
# Wang Ming wangmcas(AT)gmail.com
# 2015-01-23
#################################################

use strict;
use warnings;
use LWP::UserAgent;
use POSIX qw(strftime);
use File::Basename qw(dirname basename);
use File::Spec::Functions qw(catdir catfile);
use Cwd qw(abs_path cwd);
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

## Create working dir
my $in_gff = '';
my $out_dir = 'tmp';
my $help = '';
my $man = 0;
my $section = 100;

GetOptions(
        'input|i=s' => \$in_gff,
        'output|o=s' => \$out_dir,
        'section|s=i' => \$section,
        'help|h' => \$help,
        'man' => \$man
        ) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitval => 0, -verbose => 2) if $man;

#die 'Need -i input.gff' unless($in_gff);

## Prepare working dir : $out_dir
if(-d $out_dir){
    warn 'The output dir exists: ', $out_dir, "\n";
}else{
    mkdir $out_dir;
}

###################################
############# Start ###############
###################################
my $working_dir = cwd;

my $result_dir = catdir($working_dir, $out_dir);

###### PRE ########################
my $dash_line = '-' x 40;
my $dtime = &getTIME();
print 'Retrieve GO annotation of gene by GFF file'. "\n". $dash_line . "\n";
print $dtime. "\n\t". 'Step 1. Parsing GFF file: '. $in_gff. "\n";

### 1. Trans GFF to gene table#####
my $gff_basename = basename($in_gff);
$gff_basename =~ s/\.gff$/\.table/;
my $geneTable_file = catfile($result_dir, $gff_basename);
open GFF, ">", $geneTable_file or die $!;
my $geneTable = &GFF2GeneTable($in_gff);
print GFF $geneTable;
close GFF;

#### Checkpoint 1 ####
my @check_1 = split ("\n", $geneTable);
my $check_1_lines = scalar(@check_1) - 1; # exclude the header line
$dtime = &getTIME();
print $dtime. "\n\t". '- Find the number of genes: '. $check_1_lines. "\n";
print "\t". '- Writing geneTable to file: '. $geneTable_file. "\n";

## 2. Trans geneID to UniProt id
print "\t". 'Step 2. Retrieve UniProt id of gene by geneID/entrezgene id'. "\n";
my $UniProt_base = basename($in_gff);
$UniProt_base =~ s/\.gff$/\.UniProt.txt/;
my $UniProt_file = catfile($result_dir, $UniProt_base);
open UNIP, ">", $UniProt_file or die $!;
my @geneIDs = &GeneTable2IDs($geneTable);
my $geneIDs_length = scalar(@geneIDs);
my %UP2GI = ();

#my @query_ids = &splitArrays(' ', '1000', @geneIDs);
my @query_ids = &splitArrays(' ', $section, @geneIDs);

foreach my $query_id (@query_ids){
    my $Input_format = 'P_ENTREZGENEID';
    my $Output_format = 'ACC';
    my $Input_ids = $query_id;    
    my %res = &GeneID2UniProt($Input_format, $Output_format, $Input_ids);
    foreach my $k (keys %res){ # $k = GeneID, value = UniProt id
        $UP2GI{$res{$k}} = $k;
        print UNIP $k. "\t". $res{$k}. "\n";;
    }
}

close UNIP;

#### Checkpoint 2 ####
my @check_2 = keys %UP2GI;
my $check_2_lines = scalar(@check_2);
$dtime = &getTIME();
print $dtime. "\n\t". '- Find UniProt ids: '. $check_2_lines. ' of '. $check_1_lines. "\n";
print "\t". '- Writing UniProt results to file: '. $UniProt_file. "\n";

## 3. Retrieve GO annotation by UniProt id
print "\t". 'Step 3. Retrieve GO annotation by UniPort id'. "\n";
my $go_anno_base = basename($in_gff);
$go_anno_base =~ s/\.gff$//;
my $go_anno_file = catfile($result_dir, "$go_anno_base\.GOanno.txt");
my $go_entrez_index = catfile($result_dir, "$go_anno_base\.GO.db");
open ANNO, ">", $go_anno_file or die $!;
open LIST, ">", $go_entrez_index or die $!;
my @UniP_ids = keys %UP2GI;
my $UniP_length = scalar(@UniP_ids);
my %check_3 = ();

#my @query_UniProts = &splitArrays(',', '1000', @UniP_ids);
my @query_UniProts = &splitArrays(',', $section, @UniP_ids);
my %check_Num = ();

foreach my $query_UniProt (@query_UniProts){
    my $GOdb = &UniProt2GO($query_UniProt);
    my $smpGO = &smpGOanno($GOdb);
    print ANNO $GOdb;
    print LIST $smpGO;
    push @{$check_3{'anno'}}, (split "\n", $GOdb);
    push @{$check_3{'smp'}}, (split "\n", $smpGO);
    my @lines = split (/\n/, $smpGO);
    foreach my $l(@lines){
        my ($u, $g) = split /\t/,$l;
        $check_Num{'UniProt'}->{$u} = 1;
        $check_Num{'GO'}->{$g} = 1;
    }
}
my $Num_UniProt = scalar(keys %{$check_Num{'UniProt'}});
my $Num_GO = scalar(keys %{$check_Num{'GO'}});
close ANNO;
close LIST;

#### Checkpoint 3 ####
my $check_3_anno = scalar(@{$check_3{'anno'}});
my $check_3_smp = scalar(@{$check_3{'smp'}});
$dtime = &getTIME();
print $dtime. "\n\t". '- Retrieve the number of GO annotations: '. $check_3_anno. "\n";
print "\t". '- Ouput entrezgene+go_id lines: '. $check_3_smp. "\n";
print "\t". '- Retrieved '. $Num_UniProt. ' UniProt ids, with '. $Num_GO. ' of unique GO annotations.'. "\n";
print "\t". '- Writing full go annotations to file: '. $go_anno_file. "\n";
print "\t". '- Writing entrezgene + GO_id to file: '. $go_entrez_index. "\n";
print $dtime. "\n\t". 'Finished.'. "\n";
print $dash_line. "\n";

#### Clear temp files ####
unlink glob "annotation-*.tsv";

######## SUBROUTINE #########
sub GFF2GeneTable{
    my $in_gff = shift(@_);
    open F, $in_gff or die $!;
    my $header = join "\t", ('GeneID', 'seqname', 'start', 'end', 'strand', 'GeneName', 'Locus');
    my $genetable .= $header . "\n";
    while(<F>){
        chomp;
        next unless(/\tgene\t/);
        my ($seqname, $start, $end, $strand, $description) = (split(/\t/, $_))[0,3,4,6,8];
        my ($GeneName, $GeneID, $Locus) = qw(- - -);
        ($GeneName) = ($description =~ /locus_tag=(\w+)/);
        if($description =~ /Name.*gbkey.*locus_tag/){ # gff 1.20
            ($GeneName, $GeneID, $Locus) = ($_ =~ /Name=(.*)\;Dbxref=GeneID\:(\d+)\;.*locus_tag=(\w+)/);
            ($GeneName) = ($GeneName =~ /note=(\w+\-\w+)/)?$1:$GeneName;
        }else{ # gff 1.14
            ($Locus, $GeneID) = ($description =~ /locus_tag=(\w+)\;.*GeneID:(\d+)/);
            if($description =~ /ID=\w+.\d\:(\w+)\;locus_tag/){
                $GeneName = $1;
            }else{
                $GeneName = '-';
            }
        }
        $GeneName = (split /\;/, $GeneName)[0] if($GeneName =~ /\;/);
        my $line = join "\t", ($GeneID, $seqname, $start, $end, $strand, $GeneName, $Locus);
        $genetable .= $line . "\n";
    }
    return $genetable;
}

### GeneID2UniProt
sub GeneID2UniProt{
    my ($from, $to, $input_ids) = @_;
    my $base = 'http://www.uniprot.org';
    my $tool = 'mapping';
    my $params = {
        from => $from,
        to => $to,
        format => 'tab',
        query => $input_ids
    };
    #Please set your email address here to help us debug in case of problems.
    my $contact = 'wangmcas@gmail.com';
    my $agent = LWP::UserAgent->new(agent => "libwww-perl $contact");
    push @{$agent->requests_redirectable}, 'POST';

    my $response = $agent->post("$base/$tool/", $params);
    while (my $wait = $response->header('Retry-After')) {
        print STDERR "Waiting ($wait)...\n";
        sleep $wait;
        $response = $agent->get($response->base);
    }
    my %result;
    if($response->is_success){
        my @lines = split (/\n/, $response->content);
        foreach my $line (@lines){
            next unless( my ($f, $t) = split /\s+/, $line);
            $result{$f} = $t if $f ne 'From';
        }
    }else{
        die 'Failed, got ' . $response->status_line . ' for ' .
            $response->request->uri . "\n";
    }
    return %result;
}

### GeneTable to IDs ###
sub GeneTable2IDs{
    my $table = shift(@_);
    my @lines = split (/\n/,$table);
    my @IDs = ();
    foreach my $l (@lines){
        next if($l =~ /^GeneID/i);
        my $id = (split /\s+/, $l)[0];
        push @IDs, $id;
    }
    return @IDs;
}

### UniProt to GO annotation ###
sub UniProt2GO{
    my $querys = shift(@_);
    my $ua = LWP::UserAgent->new;
    my $url = 'http://www.ebi.ac.uk/QuickGO/GAnnotation?protein=';
    my $fmt = '&format=tsv&col=proteinID,proteinSymbol,evidence,goID,goName,aspect,ref,with,from';
    my $get = $url . $querys . $fmt;
    my $req = HTTP::Request->new(
            GET => $get
            );
    my $annotation_tsv = 'annotation-'. $$ . '.tsv';
    my $res = $ua->request($req, $annotation_tsv);
    open (FILE, $annotation_tsv);
    my $head = <FILE>;
    my $output = '';
    while(<FILE>){
        chomp;
        my ($proteinID, $proteinSymbol, $evidence, $goID, $goName, $aspect, $ref, $with, $from) = split(/\t/);
        $output .= "$proteinID => $goID ($aspect: $goName) $evidence $ref $with $from\n";
    }
    close FILE;
    return $output;
}

### simplify GO annotation output ###
sub smpGOanno{
    my $inGO = shift(@_);
    my @lines = split ("\n", $inGO);
    my $smpGO = '';
    foreach my $l (@lines){
        my ($uniprot, $go_id) = (split /\s+/, $l)[0,2];
        my $entrezgene = ($UP2GI{$uniprot})?$UP2GI{$uniprot}:'-';
#        my $entrezgene = $UP2GI{$uniprot};
        $smpGO .= $entrezgene. "\t". $go_id. "\n";
    }
    return $smpGO;
}

### split input Array ###
sub splitArrays{ # ',' 400, @arrays
    my $sep = shift(@_);
    my $num = shift(@_);
    my @in_arrays = @_;
    my $in_length = scalar(@in_arrays);
    my @out_arrays = ();
    
    for(my $i=0; $i<=$in_length; $i+=$num){
        my @tags = ();
        for(my $t=$i; $t<$i+$num; $t++){
            last if($t >= $in_length);
            push @tags, $in_arrays[$t];
        }
        push @out_arrays, join ($sep, @tags);
    }
    return @out_arrays;
}

### get localtime ###
sub getTIME{
    my $datestring = strftime "%Y-%m-%d %H:%M:%S", localtime;
    return $datestring;
}

__END__

=head1 NAME

c<GFF2GOannotation.pl> - Retrieve GO annotation of genes in GFF file.

=head1 SYNOPSIS

perl GFF2GOannotation.pl -i NC_000962.gff -o Rv_GOmap -s 400

=head1 OPTIONS

=over 8

=item B<-i> <-input>

The GFF file in GFF3 format. (bacteria, one chromosome).

=item B<-o> <-output>

The output directory for all the results. default [tmp]

=item B<-s> <-section>

The number of queries for UniPort Web Service and quickGO Web Services.
default [100]

=item B<-h>

Print this help.

=item B<-man>

Show more details.

=back

=head1 DESCRIPTION

This scirpt consist of 3 steps of transformation to get gene GO annotation from quickGO.  

1. Parse the GeneID (entrezgene) id of each gene/CDS from GFF file.  

2. Retrieve the UniProt/Swiss-Prot ID of each gene by (entrezgene) id from 
UinProt/Swiss-Prot Web Services.   
(http://www.uniprot.org/help/programmatic_access)  

3. Retrieve the GO annotations for each gene by (UniProt) id from 
quickGO Web Services.     
(http://www.ebi.ac.uk/QuickGO/WebServices.html)  

Produce a 2-column table: "entrezgene" and "go_accession".

Output:

1. *.GO.db: A 2-column file with "entrezgene id" and "GO annotation id".

2. *.table: The gene table of input GFF file.

3. *.UniProt.txt: This file contains "entrezgene id" and "UniProt id" for each gene.

4. *.GOanno.txt: This file contains "entrezgene id", "UniProt di" and "GO annotation", "GO description".

=head1 AUTHOR

Wang Ming wangmcas<at>gmail.com

Last Update: Feb 2, 2015

=cut

