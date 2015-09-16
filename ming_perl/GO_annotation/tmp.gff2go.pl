#!/usr/bin/env perl

#################################################
# Find GO annotation of genes (from GFF file)
#
# Wang Ming wangmcas(AT)gmail.com
# 2015-01-23
#################################################

use strict;
use warnings;
use LWP::UserAgent;
use File::Basename qw(dirname basename);
use File::Spec::Functions qw(catdir catfile);
use File::Path qw(make_path remove_tree);
use Cwd qw(abs_path cwd);
use POSIX qw(strftime);
use Getopt::Std;
use Data::Dumper;

gff2go();
exit(1);

sub gff2go {
    my %opts = (n => 400);
    getopts("o:i:n:u:h", \%opts);
    usage() if(defined $opts{h});
#    usage_full() if(defined $opts{});
    die("[-i] GFF file required\n") if(! defined $opts{i});
    die("[-o] output directory required\n") if(! defined $opts{o});
    print 'Retrieve GO annotations of the genes from GFF file' . "\n" . '-' x 40 . "\n";
    ### prepare working dir
    my $wk_dir  = cwd;
    my $out_dir = catdir($wk_dir, $opts{o});
    make_path($out_dir) if( ! -d $out_dir);

    print  show_date() . ' - Parsing genes from GFF file:' . basename($opts{i}) . "\n";
    ### parsing GFF file 
    my %tab = %{ gff2table($opts{i}) };
    my $tab_file = catfile($out_dir, basename($opts{i}));
    $tab_file =~ s/.gff$/.table/;
    open my $fh_tab, "> $tab_file" or die "Cannot open $tab_file, $!\n";
    print $fh_tab join("\t", 'GeneID', 'seqname', 'start', 'end', 'strand', 'GeneName', 'Locus') . "\n";
    print $fh_tab join("", values %tab);
    close $fh_tab;
    my @entrez_ids = keys %tab;
    
    # checkpoint: number of hits 
    die("[-i $opts{i}] genes not found in GFF file\n") if((keys %tab) < 2);
    print show_date() . ' - Save gene table to file: ' . basename($tab_file) . ' [' . scalar(@entrez_ids) . "\]\n";
    ### start
    my %e2u = ();
    my $uniprot_list = catfile($out_dir, basename($opts{i}));
    $uniprot_list =~ s/\.gff/.Uniprot.txt/;
    if(defined $opts{u}) {
        ### Download proteomes annotation of input genome
        ### ouput as Tab-separated format (7-column)
        ### [Entry Entry_name Status Protin_names Gene_names Organism Length]
        ### multiple genes can be found in [Gene_names] field.
        %e2u = %{ addNotes($tab_file, $opts{u}, $uniprot_list) }; 
    }else{
        ### Using UniProt API to retrieve ids
        ### WARN: may miss some queries for unknown reason.
        my @entrez_id_chunks = chunk_ids(\@entrez_ids, $opts{n}, " ");
        for my $e (@entrez_id_chunks) {
            my %sub = %{ entrez2uniprot($e) };
            %e2u = (%e2u, %sub);
        }
        open my $fh_up, "> $uniprot_list" or die "Cannot open $uniprot_list, $!\n";
        for my $u (sort keys %e2u) {
            for my $e (@{$e2u{$u}}) {
                print $fh_up, join("\t", $e, $u) . "\n";
            }
        }
        close $fh_up;
    }
    # checkpoint: number of UniProt hits
    die("[No UniProt ids found] \n") if((keys %e2u) < 1);
    print show_date() . ' - Save UniProt ids to file: ' . basename($uniprot_list) . ' ['. scalar(keys %e2u) . "\]\n";

    ### UniProt to GO 
    my $go_list = catfile($out_dir, basename($opts{i}));
    $go_list =~ s/\.gff/.GO.anno.txt/;
    open my $fh_st, "> $go_list" or die "Cannot open $go_list, $!\n";
    my @uniprot_ids = keys %e2u;
    my @ui_chunks = chunk_ids(\@uniprot_ids, $opts{n}, ",");
    my $u2g = '';
    for my $ui (@ui_chunks) {
        my $an = uniprot2go($ui, $out_dir);
        next if($an eq ''); # do not return results
        print $fh_st $an;
        $u2g .= $an;
    }
    close $fh_st;

    #checkpoint: number of GO hits
    my @go_lines = split /\n/, $u2g;
    die("[Not found GO annotations\n]") if(@go_lines < 1);
    print show_date() . ' - Save GO annotation to file: ' . basename($go_list) . ' [' . scalar(@go_lines) . "\]\n";
    
    my %report = ();
    ### create GO.db
    my $ez2go = '';
    for my $g (@go_lines) {
        my ($n_un, $n_go) = (split /\s+/, $g)[0, 2];
        die("[UniProt ids not fount: $n_un\n]") if(! exists $e2u{$n_un});
        $report{un}->{$n_un} = 1;
        $report{go}->{$n_go} = 1;
        for my $e (@{$e2u{$n_un}}) {
            $report{ez}->{$e} = 1;
            $report{db} ++;
#            $ez2go .= $e . "\t" . $n_go . "\n";
            $ez2go .= $n_go . "\t" . $e . "\n"; # GO+GeneID, order changed in clusterProfiler v2.2.3
        }
    }
    my $go_db = catfile($out_dir, basename($opts{i}));
    $go_db   =~ s/\.gff/.GO.db/;
    open my $fh_db, "> $go_db" or die "Cannot open $go_db, $!\n";
    print $fh_db $ez2go;
    close $fh_db;
    print show_date() . ' - Save GO db to file: ' . basename($go_db) . ' [' . $report{db} . "\]\n";
    print show_date() . ' - Finish' . "\n" . '-' x 40 . "\n";
    ### summary
    print 'Summary' . "\n";
    printf "%-15s: %-8s\n", 'Total genes', scalar(@entrez_ids);
    printf "%-15s: %-8s\n", 'Genes in db', scalar(keys %{$report{ez}});
    printf "%-15s: %-8s\n", 'UniPort in db', scalar(keys %{$report{un}});
    printf "%-15s: %-8s\n", 'GO in db', scalar(keys %{$report{go}});
}

### split array by specific length
sub chunk_ids {
    my @ids = @{$_[0]};
    my $num = $_[1];
    my $sep = $_[2]; # " ", ",", ";"
    my @out = ();
    for(my $i = 0; $i <= $#ids; $i += $num) {
        my @sub = ();
        for(my $k = $i; $k <= $i + $num; $k ++) {
            last if($k > $#ids);
            push @sub, $ids[$k];
        }
        push @out, join($sep, @sub);
    }
    return(@out);
}

### Read standard NCBI GFF file
sub gff2table {
    my $in = $_[0];
    my %tb = ();
    open my $fh_in, "< $in" or die "Cannot open $in, $!\n";
    while(<$fh_in>) {
        chomp;
        next if(/^\#|^\s*$/);
        next unless(/\tgene\t/);
        my ($id, $start, $end, $strand, $des) = (split /\t/, $_)[0, 3, 4, 6, 8];
        my $name    = '-';
        my $gene_id = '-';
        my $locus   = '-';
        if(($name)    = $des =~ /Name=(\w+)/) {
        }elsif($name  = $des =~ /ID=(\w+)/) {}
        if(($gene_id) = $des =~ /GeneID\:(\d+)/) {}
        if(($locus)   = $des =~ /locus_tag=(\w+)/) {}
        $name = (split /\;/, $name)[0] if($name =~ /\;/);
        $tb{$gene_id} = join("\t", $gene_id, $id, $start, $end, $strand, $name, $locus) . "\n";
    }
    return (\%tb);
} 

sub addNotes {
    my $f1   = $_[0];
    my $f2   = $_[1];
    my $out  = $_[2];
    ###
    my %h = ();
    open my $fh_f2, "< $f2" or die "Cannot open $f2, $!\n";
    while(<$fh_f2>) {
        chomp;
        my @tabs  = split /\t/;
        die("[-u $f2] is not correct, see: http://www.uniprot.org/proteomes \n") if(@tabs != 7);
        my @genes = split /\s+/, $tabs[4];
        for my $g (@genes) {
            push @{$h{$g}}, $tabs[0];
        }
    }
    close $fh_f2;
    ###
    my %uns = ();
    open my $fh_f1, "< $f1" or die "Cannot open $f1, $!\n";
    open my $fh_st, "> $out" or die "Cannot write to $out, $!\n";
    while(<$fh_f1>) {
        chomp;
        my @tabs = split /\t/;
        my $note = '-';
        if(exists $h{$tabs[6]}) {
            my @notes = ();
            for my $g (@{$h{$tabs[6]}}) {
                push @{$uns{$g}}, $tabs[0];
                push @notes, $g;
            }
            $note = join(',', @notes);
        }
        print $fh_st join("\t", @tabs, $note) . "\n";
    }
    close $fh_f1;
    close $fh_st;
    ###
    return(\%uns);
}

### EntrzeID to UniProt id
### from P_ENTREZGENEID to ACC
### find details at: http://www.uniprot.org/help/programmatic_access#batch_retrieval_of_entries
sub entrez2uniprot {
    my $id  = $_[0]; # entrez_id
    ### Using UniProt REST
    my $base = 'http://www.uniprot.org';
    my $tool = 'mapping';
    my $params = {
        from => 'P_ENTREZGENEID',
        to   => 'ACC',
        format => 'tab',
        query  => "$id"  # eg: 14515852 14515879 14515880
    };
    my $contact = 'wangmcas@gmail.com';
    my $agent = LWP::UserAgent->new(agent => "libwww-perl $contact");
    push @{$agent->requests_redirectable}, 'POST';
    my $response = $agent->post("$base/$tool/", $params);
    while (my $wait = $response->header('Retry-After')) {
        print STDERR "Waiting ($wait)...\n";
        sleep $wait;
        $response = $agent->get($response->base);
    }
    my %hit = ();
    if($response->is_success){
        my @lines = split /\n/, $response->content;
        for my $n (@lines) {
            next if($n =~ /From\s+To/);
            my ($ez_id, $un_id) = split /\s+/, $n;
            push @{$hit{$un_id}}, $ez_id;
        }
    }else{
        die('Failed, got ' . $response->status_line . ' for ' . 
            $response->request->uri . "\n");
    }
    return (\%hit);
}

### UniProt id to GO annotation
### Find more details at: http://www.ebi.ac.uk/QuickGO/WebServices.html
### search all GO using uniprot ids
sub uniprot2go {
    my $query_id  = $_[0];  # multiple ids need to be connected by ',': P12345,Q4VCS5
    my $outdir    = $_[1];
    my $ua  = LWP::UserAgent->new;
    my $url = 'http://www.ebi.ac.uk/QuickGO/GAnnotation?protein=';
    my $fmt = '&format=tsv&limit=-1&col=proteinID,proteinSymbol,evidence,goID,goName,aspect,ref,with,from';
    my $get = $url . $query_id . $fmt;
    my $req = HTTP::Request->new( GET => $get );
    my $anno_tsv = catfile($outdir, 'anno-' . $$ . '.tsv');
    my $res = $ua->request($req, $anno_tsv);
    open my $fh_tsv, "< $anno_tsv" or die "Cannot open $anno_tsv, $!\n";
    my $head = <$fh_tsv>;
    my $out;
    while(<$fh_tsv>){
        chomp;
        my ($proteinID, $proteinSymbol, $evidence, $goID, $goName, $aspect, $ref, $with, $from) = split//, '-' x 9;
        ($proteinID, $proteinSymbol, $evidence, $goID, $goName, $aspect, $ref, $with, $from) = split(/\t/);
        $out .= "$proteinID => $goID ($aspect: $goName) $evidence $ref $with $from\n";
    }
    close $fh_tsv;
    return $out;
}


sub show_date {
    my $dt = strftime "[%Y-%m-%d %H:%M:%S]", localtime;
    return $dt;
}

#################################################
sub usage {
    die("
Usage: gff2go.pl [options] <in.gff>

Options: -o <STR>       : redirect the results to a directory
         -i <STR>       : GFF file
         -u <STR>       : Input the UniProt annotations 
                          (recommended, download UniProt ids online:)
         -n <INT>       : number of items submit to UniProt at once [400]

Example:
gff2go.pl -o out -i NC_000962.gff -u uniport-ids.tab
\n");
}

sub usage_full {
    die('
Usage: gff2go.pl [options] <in.gff>

Options: -o <STR>       : redirect the results to a directory
         -i <STR>       : GFF file
         -u <STR>       : Input the UniProt annotations 
                          (recommended, download UniProt ids online:)
         -n <INT>       : number of items submit to UniProt at once [400]

Example:
gff2go.pl -o out -i NC_000962.gff -u uniport-ids.tab

Description:
This scirpt will retrieve the GO annotation from quickGO.  
1. Parse the GeneID (entrezgene) id of each gene/CDS from GFF file.  
2. Retrieve the UniProt/Swiss-Prot ID of each gene by (entrezgene) id from 
   UinProt/Swiss-Prot Web Services.   
   (http://www.uniprot.org/help/programmatic_access)  
3. Retrieve the GO annotations for each gene by (UniProt) id from 
   quickGO Web Services.     
   (http://www.ebi.ac.uk/QuickGO/WebServices.html)  

Output:
1. *.GO.db: A 2-column file with "entrezgene id" and "GO annotation id".
2. *.table: The gene table of input GFF file.
3. *.UniProt.txt: This file contains "entrezgene id" and "UniProt id" for each gene.
4. *.GOanno.txt: This file contains "entrezgene id", "UniProt di" and "GO annotation", "GO description".

Wang Ming wangmcas(AT)gmail.com
Created: Feb 2, 2015
Modified: July 17, 2015
');
}
