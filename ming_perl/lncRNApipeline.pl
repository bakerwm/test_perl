#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Temp;
use File::Basename;
use Cwd qw (abs_path cwd);
use Term::ANSIColor qw(:constants);
Getopt::Long::Configure('prefix=--');

my ($LASTCHANGEDBY) = q$LastChangedBy: konganti $ =~ m/.+?\:(.+)/;
my ($LASTCHANGEDDATE) = q$LastChangedDate: 2015-07-08 11:45:27 -0500 (Fri, 08 May 2015)  $ =~ m/.+?\:(.+)/;
my ($VERSION) = q$LastChangedRevision: 0.7.9 $ =~ m/.+?\:\s*(.*)\s*.*/;
my $AUTHORFULLNAME = 'Kranti Konganti';

my ($help, $quiet, $setup, $get_uq_sc_opts,
    $fetch_sc_opts, $categorize_sc_opts,
    $cuffcmp_opts, $deps, $out, $start_cpc, 
    $start_rnafold, $rm_int_plots, $start_infernal,
    $num_cpu, $local_lib, $inf_cov, $skip_cpc_core,
    $skip_rnafold_core, $skip_cmscan_core,
    $lncRNApipe_ver_info, $no_update_check, 
    $skip_get_uq, $lncRNApipe_succ, $num_cpu_by_2,
    $skip_cpc, $debug);

my $is_valid_option = GetOptions('help:s'                           => \$help,
                                 'quiet'                            => \$quiet,
				 'setup'                            => \$setup,
				 'run=s'                            => \$out, 
				 'cuffcompare|cuffcmp=s'            => \$cuffcmp_opts,
				 'cat-ncRNAs=s'                     => \$categorize_sc_opts,
				 'fetch-seq:s'                      => \$fetch_sc_opts,
				 'get-uq-feat=s'                    => \$get_uq_sc_opts,
				 'cpc'                              => \$start_cpc,
				 'rnafold:s'                        => \$start_rnafold,
				 'rm-int-plots'                     => \$rm_int_plots,
				 'infernal:s'                       => \$start_infernal,
				 'cpu=i'                            => \$num_cpu,
				 'skip-blastall-core'               => \$skip_cpc_core,
				 'skip-rnafold-core'                => \$skip_rnafold_core,
				 'skip-cmscan-core'                 => \$skip_cmscan_core,
				 'coverage-infernal|cov-inf=f'      => \$inf_cov,
				 'local-lib=s'                      => \$local_lib,
				 'version'                          => \$lncRNApipe_ver_info,
				 'no-update-chk|noup|d'             => \$no_update_check,
				 'skip-get|skip-get-uq'             => \$skip_get_uq,
				 'skip-cpc'                         => \$skip_cpc,
				 'send-run-report=s'                => \$debug);

# Initialize known defaults
my $CAT_SC = 'categorize_ncRNAs.pl';
my $FETCH_SC = 'fetch_seq_from_ucsc.pl';
my $GET_UQ_SC = 'get_unique_features.pl';
my $CPC_SH_SC = 'CPC.sh';
my $ANN_INF_SC = 'annotate_final_ncRNAs.sh';
my $RELPLOT_MOD_SC = 'relplot_mod.pl';
my $SEND_MAIL_SC = 'send_lncRNApipe_log.pl';
my $USER_HOME = dirname(abs_path($0));
if (defined $num_cpu) {
    $num_cpu_by_2 = sprintf("%.0f", $num_cpu / 2) - 1;
} 
else {
    $num_cpu_by_2 = 1;
}
$inf_cov = 10 if (!defined $inf_cov);
$skip_get_uq = 1 if (defined $skip_get_uq || (!defined $get_uq_sc_opts && defined $fetch_sc_opts));

if (!$is_valid_option) {
    die "\nInvalid / No Option(s) Provided:$! See $0 --help for valid options.\n\n";
}
elsif (defined $lncRNApipe_ver_info) {
    my $io = run_lncRNApipe('version');
    $io->this_script_info('lncRNApipe',
			  $VERSION,
			  $AUTHORFULLNAME,
			  $LASTCHANGEDBY,
			  $LASTCHANGEDDATE, 'gnu');
    exit 0;
}
elsif (defined $setup) {
    setup_lncRNApipe(`pwd`) 
}
elsif (defined $out || (defined $help && $help eq '')) {
    run_lncRNApipe();
}
elsif (defined $debug) {
    my $depconf_fh = read_conf('/.lncRNApipe.depconf');
    my $pl_dep_paths_fh = read_conf('/.lncRNApipe.PERLLIBS');
    my $pl_inc_string = source_env($pl_dep_paths_fh);
    get_deps($depconf_fh);
    system("$deps->{'email'} $debug");
}
elsif (defined $help && $help ne '') {
    my $depconf_fh = read_conf('/.lncRNApipe.depconf');
    my $pl_dep_paths_fh = read_conf('/.lncRNApipe.PERLLIBS');
    my $pl_inc_string = source_env($pl_dep_paths_fh);
    get_deps($depconf_fh);
    system("$deps->{'cuffcompare'} -h") if ($help =~ m/cuff|cuffcompare/i);
    system("perl $pl_inc_string $deps->{'cat'} -h") if ($help =~ m/cat|cat-ncRNAs/i);
    system("perl $pl_inc_string $deps->{'get'} -h") if ($help =~ m/get|get-uq-feat/i);
    system("perl $pl_inc_string $deps->{'fetch'} -h") if ($help =~ m/fetch|fetch-seq/i);
    system("$deps->{'RNAfold'} --detailed-help") if ($help =~ m/rnafold|rna|fold/i);
    system("$deps->{'cmscan'} -h") if ($help =~m /infernal|inf/i);
    print "\nSee http://cpc.cbi.pku.edu.cn/docs/install_guide.jsp for help documentation.\n\n" if ($help =~ m/cpc/i);
    close $depconf_fh;
    close $pl_dep_paths_fh;
}
else {

################################################################################################
#                                                                                              #
#     Check for updates first, unless disabled as some serious bugs may have been fixed        #
#     =================================================================================        #
#                                                                                              #
################################################################################################
    check_for_updates();
}

################################################################################################
#                                                                                              #
#                                     Run the pipeline                                         #
#                                     ================                                         #
#                                                                                              #
################################################################################################

sub run_lncRNApipe {
    my $shall_I_run = shift;
    my $outdir;
    my $pl_dep_paths_fh = read_conf('/.lncRNApipe.PERLLIBS');
    my $pl_inc_string = source_env($pl_dep_paths_fh);
    $shall_I_run = '' if (!defined $shall_I_run);

    # We will quit, if the setup was unsuccessful.
    require IO::Routine;
    IO::Routine->import();

    my $io = IO::Routine->new($help, $quiet);

    return $io if ($shall_I_run eq 'version');

    $io->verify_options([$out, $is_valid_option]) if (defined $help && $help eq '');

    if (!defined $out) {
	$outdir = $io->validate_create_path($out, 'do not create',
					    'lncRNApipe pipeline output directory');
    }
    else {
	$outdir = $io->validate_create_path($out, 'create',
					    'lncRNApipe pipeline output directory')
    }

    my $s_time = $io->start_timer;
    binmode(STDOUT, ":utf8");

    # Validate options here after help documentation is printed
    $io->c_time('Validating options...');
    $io->c_time("Starting \x{2632}\x{2634} lncRNApipe Pipeline...");

    my $dep_tools_fh = $io->open_file('<', $USER_HOME . '/.lncRNApipe.depconf');
    get_deps($dep_tools_fh);

    $io->exist_sys_cmd(['mkdir']);

    my $cuffcmp_dir = $outdir . 'cuffcompare';    

    # Run cuffcompare ################################################################################
    do {

	module_header($io, 'Module 1: Running cuffcompare...');
	
	$io->execute_system_command("mkdir -p $cuffcmp_dir",
				    'Making output directory for cuffcompare [ ' . $cuffcmp_dir . ' ]')
	    if (!-d $cuffcmp_dir);

	# Check if -i is provided with cuffcompare and adjust accordingly.

	(my $cuffcmp_opts_trs_in_file) = ($cuffcmp_opts =~ m/-i\s+([^\s]*)/);
	
	if ($cuffcmp_opts_trs_in_file) {
	    my $i_fh = $io->open_file('<', $cuffcmp_opts_trs_in_file);
	    my $i_files = do {local $/; <$i_fh>};
	    my $i_file_str = join(' ', split(/\n|\r/, $i_files));
	    
	    $cuffcmp_opts =~ s/-i\s+[^\s]*/$i_file_str/;
	}

	$io->execute_system_command($deps->{'cuffcompare'} .
				    " -T -o $cuffcmp_dir/lncRNApipe_cuffcmp " .
				    join('-', split(/\-/, $cuffcmp_opts)),
				    "Command call:\n" .
				    "-------------\n" . $deps->{'cuffcompare'} .
				    " -T -o $cuffcmp_dir/lncRNApipe_cuffcmp " .
					join('-', split(/\-/, $cuffcmp_opts)));

    } if (defined($cuffcmp_opts) && $cuffcmp_opts ne '');
    

    # Run categorize_ncRNAs.pl ########################################################################
    do {

	module_header($io, 'Module 2: Running categorize_ncRNAs.pl');
	
	my $is_failed = '';
	my $cat_dir = $outdir . 'categorize_ncRNAs';
	
	if (!-e "$cuffcmp_dir/lncRNApipe_cuffcmp.tracking" || !-s "$cuffcmp_dir/lncRNApipe_cuffcmp.tracking") {
	    
	    print "\nCannot find Cuffcompare tracking [ $cuffcmp_dir/lncRNApipe_cuffcmp.tracking ] file...\n";
	    $lncRNApipe_succ = pipeline_status($io, 'ERROR'); 
	}
	
	$io->execute_system_command("mkdir -p $cat_dir",
				    'Making output directory for ' .
				    $io->file_basename($deps->{'cat'}, 'suffix') .
				    ' [ ' . $cat_dir . ' ]')
	    if (!-d $cat_dir);
	
	if (defined $cuffcmp_opts && $cuffcmp_opts ne '') {

	    $lncRNApipe_succ = pipeline_status($io, 'ERROR!', "Required options not mentioned for -cat module.\nSee '" .
					       "$0 -h cat" . "' for required options.") if ($categorize_sc_opts !~ m/sample/);
	    
	    my ($tr_files, $tr_size) = get_num_tr_files($cuffcmp_opts);

	    my $annot_for_cat = '';
	    if ($categorize_sc_opts !~ m/-ann.*?\s+.+/) { 
		my $annot = get_annot($cuffcmp_opts);
		
		$io->c_time('Converting annotation file supplied with cuffcompare to gene prediction format [ GTF -> genePred ]');
		$annot_for_cat = '-annot ' . $cat_dir . '/' . $io->file_basename($annot) . '.txt';
		
		my $gtfToGenePred_failed = $is_failed = $io->execute_get_sys_cmd_output($deps->{'bin-gtfToGenePred'} . ' -genePredExt -geneNameAsName2 ' .
											$annot . ' ' . $cat_dir . '/' . $io->file_basename($annot) . '.txt',
											"Command call:\n" .
											"-------------\n" .
											$deps->{'bin-gtfToGenePred'} . ' -genePredExt -geneNameAsName2 ' .
											$annot . ' ' . $cat_dir . '/' . $io->file_basename($annot) . '.txt');
		
		$gtfToGenePred_failed =~ s/[\s\r\n]+//g;
		$io->error('gtfToGenePred quit with the following error:' . "\n\n$is_failed\nBailing out!"),
		$lncRNApipe_succ = 0
		    if ($gtfToGenePred_failed =~ m/invalid|error|not/i);
	    }
	    
	    $categorize_sc_opts .= " -cpu $num_cpu" if (defined $num_cpu);
	    
	    $is_failed = $io->execute_get_sys_cmd_output('perl' . $pl_inc_string . $deps->{'cat'} . ' ' .
							 join('-', split(/\-|\--/, $categorize_sc_opts)) .
							 " -cuffcmp $cuffcmp_dir/lncRNApipe_cuffcmp.tracking -out $cat_dir " .
							 "-bin $deps->{'bin-gtfToGenePred'} " . $annot_for_cat . ' ' .
							 $tr_files, 
							 "Command call:\n" .
							  "-------------\n" . $deps->{'cat'} . ' ' . 
							 join('-', split(/\-|\--/, $categorize_sc_opts)) . 
							 " -cuffcmp $cuffcmp_dir/lncRNApipe_cuffcmp.tracking -out $cat_dir " .
							 "-bin $deps->{'bin-gtfToGenePred'} " . $annot_for_cat . ' ' .
							 $tr_files);
	}
	else {
	    		        
	    $lncRNApipe_succ = pipeline_status($io, 'ERROR!', "Required options not mentioned for -cat module.\nSee '" .
					       "$0 -h cat" . "' for required options.") if ($categorize_sc_opts !~ m/cuff|annot|out|sample/);
	    
	    $deps->{'cat'} .= " -cpu $num_cpu" if (defined $num_cpu);
	    
	    $is_failed = $io->execute_get_sys_cmd_output('perl' . $pl_inc_string . $deps->{'cat'} . ' ' .
							 "-cuffcmp $cuffcmp_dir/lncRNApipe_cuffcmp.tracking -out $cat_dir " .
							 "-bin $deps->{'bin-gtfToGenePred'} " .
							 join('-', split(/\-|\--/, $categorize_sc_opts)),
							 "Command call:\n" .
							 "-------------\n" . $deps->{'cat'} . ' ' .
							 "-cuffcmp $cuffcmp_dir/lncRNApipe_cuffcmp.tracking -out $cat_dir " .
							 "-bin $deps->{'bin-gtfToGenePred'} " .
							 join('-', split(/\-|\--/, $categorize_sc_opts)));
	}

	$lncRNApipe_succ = pipeline_status($io, $is_failed);
	
    } if (defined($categorize_sc_opts) && $categorize_sc_opts ne '');

    # Run get_unique_features.pl ######################################################################
    do {

	my $cpu = do_parallel($num_cpu) if (defined $num_cpu);
	$lncRNApipe_succ = get_pipeline_status_code($cpu) if (defined $num_cpu);
	module_header($io, 'Module 3: Running get_unique_features.pl');
	
	my $is_failed = '';
        my $cat_dir = $outdir . 'get_unique_features';
        
	$io->execute_system_command("mkdir -p $cat_dir",
                                    'Making output directory for ' .
				    $io->file_basename($deps->{'get'}, 'suffix') .
				    ' [ ' . $cat_dir . ' ]')
            if (!-d $cat_dir);
	
	opendir (lncRNApipe_categorize_mod, $outdir . 'categorize_ncRNAs') ||
	    $io->error("Cannot open directory $outdir" .
		       'categorize_ncRNAs to read lncRNApipe class files.');
	    
	push my @lncRNApipe_class_files, grep {/putative\.class\.lncRNAs\.gtf$/} readdir lncRNApipe_categorize_mod;
	
	if ($#lncRNApipe_class_files <  0) {
	    $lncRNApipe_succ = pipeline_status($io, 'ERROR', "Cannot find output in $cat_dir ...");
	}

        $lncRNApipe_succ = pipeline_status($io, 'ERROR!', "Required options not mentioned for -get module.\nSee '" .
					   "$0 -h get" . "' for required options.") if ($get_uq_sc_opts !~ m/sf|cf/);
	
	foreach my $lncRNApipe_class_file (@lncRNApipe_class_files) {
	    	    
	    $cpu->start and next if (defined $num_cpu);
	    
	    $lncRNApipe_class_file = $outdir . 'categorize_ncRNAs/' . $lncRNApipe_class_file;
	    my $unique_ncRNAs = $cat_dir . '/' . $io->file_basename($lncRNApipe_class_file) . '.unique.gtf';
	    unlink $unique_ncRNAs if (-e $unique_ncRNAs);

	    my $sff_opt = '';
	    if ($get_uq_sc_opts !~ m/-sff/) {
		$io->c_time("Assuming supplied known ncRNAs to be in BED format...\n" .
			    'You can change this option with -sff. Ex: -sff gtf');
		$sff_opt = '-sff bed';
	    }
		
	    $is_failed = $io->execute_get_sys_cmd_output('perl' . $pl_inc_string . $deps->{'get'} . ' ' .
							 '-q -cf ' . $lncRNApipe_class_file . " -cff gtf $sff_opt ".
							 join('-', split(/\-|\--/, $get_uq_sc_opts)) .
							 ' -unique -out ' .
							 $unique_ncRNAs,
							 "Command call:\n" .
							 "-------------\n" . $deps->{'get'} . ' ' .
							 '-q -cf ' . $lncRNApipe_class_file . " -cff gtf $sff_opt ".
							 join('-', split(/\-|\--/, $get_uq_sc_opts)) .
							 ' -unique -out ' .
							 $unique_ncRNAs);
	    
	    $lncRNApipe_succ = pipeline_status($io, $is_failed);
	    if (-e $unique_ncRNAs && !-s $unique_ncRNAs) {
		$lncRNApipe_succ = pipeline_status($io, 'ERROR', "No unique features found. In fact, this is not an error.\n" .
						   'It just means that any putative lncRNAs that have been extracted from -cat module are all known ncRNAs.');
		$io->warning("No unique features found. In fact, this is not an error.\n" .
			     'It just means that any putative ncRNAs that have been extracted from -cat module are all known ncRNAs.');
		$io->warning('Proceeding with files from "categorize_ncRNAs" module.', 'INFO!!');
		$skip_get_uq = 1;
	    }
		
	    disp_ncRNA_counts($io, $unique_ncRNAs, 'Total number of putative unique ncRNAs');
	    $cpu->finish(0, {'lncRNApipe_succ' => $lncRNApipe_succ}) if (defined $num_cpu);
	}
	
	$cpu->wait_all_children if (defined $num_cpu);
	close lncRNApipe_categorize_mod;
		
    } if (defined($get_uq_sc_opts) && $get_uq_sc_opts ne '');

    # Find out if local FASTA mentioned with cuffcompare before running fetch_seq_from_ucsc.pl ########
    my $local_ref_seq = '';
    
    # Use local reference FASTA instead of fetching, if provided with cuffcompare
    # This needs to be outside do {}, to guarantee parsing, even if -fetch-seq does not have any options.
    ($local_ref_seq) = ($cuffcmp_opts =~ m/-s\s+([^\s]*)/) if ($cuffcmp_opts);	
    if ($local_ref_seq && -s $local_ref_seq && -e $local_ref_seq) {
	$io->warning('Using local FASTA reference [ ' . $io->file_basename($local_ref_seq, 'suffix') . ' ] to fetch sequences to maintain consistency.' .
			 "\n** Requires Bio::SeqIO module to be installed and available **", 'INFO!!');
	$fetch_sc_opts = "-local $local_ref_seq";
    }
    elsif( ($local_ref_seq eq '' || !-s $local_ref_seq || !-e $local_ref_seq) && !defined $fetch_sc_opts) {
	$io->error("We need to know where to fetch FASTA sequences from [ -local or -db or -ncbi ] ??\n\nSee $0 --h fetch for options.");
    }
    
    # Run fetch_seq_from_ucsc.pl ######################################################################
    do {
		
	my $cpu = do_parallel($num_cpu) if (defined $num_cpu);
	$lncRNApipe_succ = get_pipeline_status_code($cpu) if (defined $num_cpu);
	my $uq_tr_feat_dir = 'get_unique_features';
	my $tr_feat_here_suffix = 'putative\.class\.lncRNAs\.unique\.gtf$';

	module_header($io, 'Module 4: Running fetch_seq_from_ucsc.pl');
	
	my $is_failed = '';
        my $cat_dir = $outdir . 'fetch_seq_from_ucsc';
        
	$io->execute_system_command("mkdir -p $cat_dir",
                                    'Making output directory for ' . $io->file_basename($deps->{'fetch'}, 'suffix') . ' [ ' . $cat_dir . ' ]')
            if (!-d $cat_dir);

	if (defined $skip_get_uq) {
	    $uq_tr_feat_dir = 'categorize_ncRNAs';
	    $tr_feat_here_suffix = 'putative\.class\.lncRNAs\.gtf$';
	}
	
	opendir (lncRNApipe_get_unique_mod, $outdir . $uq_tr_feat_dir) ||
	    $io->error("Cannot open directory $outdir" . $uq_tr_feat_dir .
		       ' to read unique ncRNA features.');
	
	push my @lncRNApipe_unique_files, grep {/$tr_feat_here_suffix$/} readdir lncRNApipe_get_unique_mod;

	if ($#lncRNApipe_unique_files <  0) {
            $lncRNApipe_succ = pipeline_status($io, 'ERROR', "Cannot find output in $cat_dir ...");
        }
	
	foreach my $lncRNApipe_unique_file (@lncRNApipe_unique_files) {
	    
	    $cpu->start and next if (defined $num_cpu);

	    $lncRNApipe_unique_file = $outdir . $uq_tr_feat_dir . '/' . $lncRNApipe_unique_file;
	    
	    $is_failed = $io->execute_get_sys_cmd_output('perl' . $pl_inc_string . $deps->{'fetch'} . ' ' .
							 join('-', split(/\-|\--/, $fetch_sc_opts)) .
							 ' -tmap ' . $lncRNApipe_unique_file .
							 " -out $cat_dir -ff gtf -id 'transcript_id\\s+\\\"(.+?)\\\"' -skip '\\ttranscript\\t'",
							 "Command call:\n" .
							 "-------------\n" . $deps->{'fetch'} . ' ' .
							 join('-', split(/\-|\--/, $fetch_sc_opts)) .
							 ' -tmap ' . $lncRNApipe_unique_file .
							 " -out $cat_dir -ff gtf -id 'transcript_id\\s+\\\"(.+?)\\\"' -skip '\\ttranscript\\t'");
	    $lncRNApipe_succ = pipeline_status($io, $is_failed);
	    $cpu->finish(0, {'lncRNApipe_succ' => $lncRNApipe_succ}) if (defined $num_cpu);
	}
	
	$cpu->wait_all_children if (defined $num_cpu);
	close lncRNApipe_get_unique_mod;
	
    } if (defined $fetch_sc_opts);
    
    # Finally Run CPC, RNAfold and Infernal ##############################################################
    do {
	
	my $cpu = do_parallel($num_cpu) if (defined $num_cpu);
	$lncRNApipe_succ = get_pipeline_status_code($cpu) if (defined $num_cpu);

	opendir (lncRNApipe_fetch_seq_mod, $outdir . 'fetch_seq_from_ucsc') ||
	    $io->error("Cannot open directory $outdir" .
		       'fetch_seq_from_ucsc to read unique ncRNA features.');

	push my @lncRNApipe_unique_seq_files, grep {/putative\.class\.lncRNAs.*?\.fa$/} readdir lncRNApipe_fetch_seq_mod;
	my $blastall_path = sub {
	    my @file_parts = $io->file_basename(shift, 'all');
	    chop $file_parts[1] if ($file_parts[1] =~ m/\/$/);
	    return $file_parts[1]
	};
	
	$io->exist_sys_cmd(['cut']);
	
	my $lncRNApipe_final = $outdir . 'lncRNApipe.final';
	my $cat_dir = $outdir . 'CPC';
	
	$io->execute_system_command("mkdir -p $lncRNApipe_final",
				    'Making output directory to store final list of ncRNAs [ ' . $io->file_basename($lncRNApipe_final) . '  ]')
	    if (!-d $lncRNApipe_final);
	
	foreach my $lncRNApipe_unique_seq_file (@lncRNApipe_unique_seq_files) {

	    $cpu->start and next if (defined $num_cpu);
	    $num_cpu_by_2 = sprintf("%.0f", $num_cpu / scalar(@lncRNApipe_unique_seq_files)) - 1 if (defined $num_cpu);
	    $num_cpu_by_2 = 1 if ($num_cpu_by_2 == 0);

	    $lncRNApipe_unique_seq_file = $outdir . 'fetch_seq_from_ucsc/' . $lncRNApipe_unique_seq_file;
	    my $lncRNApipe_unique_seq_noncoding_file = $outdir . 'fetch_seq_from_ucsc/' . 
		$io->file_basename($lncRNApipe_unique_seq_file) . '.noncoding.FASTA';
	    my $CPC_out_file = $cat_dir . '/' . $io->file_basename($lncRNApipe_unique_seq_file) . '.CPC.predict.txt';
	    my $cpc_work_dir = $cat_dir . '/' . $io->file_basename($lncRNApipe_unique_seq_file) . '.CPC.out' ;
	    my $lncRNApipe_final_trs = $lncRNApipe_final . '/' . $io->file_basename($lncRNApipe_unique_seq_file) . '.final.gtf';
	    my $lncRNApipe_allInfernal_trs = $lncRNApipe_final . '/' . $io->file_basename($lncRNApipe_unique_seq_file) . '.allInfernalHits.gtf';
	    my $RNAfold_dir = $lncRNApipe_final . '/' . $io->file_basename($lncRNApipe_unique_seq_file) . '.pre.RNAfold';
	    my $RNAfold_dir_final = $lncRNApipe_final . '/' . $io->file_basename($lncRNApipe_unique_seq_file) . '.final.RNAfold';
	    my $RNAfold_mfe = $RNAfold_dir . '/' . $io->file_basename($lncRNApipe_unique_seq_file) . '.mfe';
	    my $unique_gtf = $outdir . 'get_unique_features/' . $io->file_basename($lncRNApipe_unique_seq_file) . '.gtf';
	    my $pre_infernal = $lncRNApipe_final . '/' . $io->file_basename($lncRNApipe_unique_seq_file) . '.pre.infernal';
	    my $final_infernal = $lncRNApipe_final . '/' . $io->file_basename($lncRNApipe_unique_seq_file) . '.final.infernal';

	    if ($skip_get_uq || !-e $unique_gtf) {
		$unique_gtf =  $outdir . 'categorize_ncRNAs/' . $io->file_basename($lncRNApipe_unique_seq_file) . '.gtf';
	    }
	    
    	    # Run CPC ####################################################################################

	    do {
		module_header($io, 'Module 5: Running CPC.sh' . ' [ on ' .  $io->file_basename($lncRNApipe_unique_seq_file, 'suffix') . ' ] ...');
	
		unlink $lncRNApipe_final_trs if (-e $lncRNApipe_final_trs);
		
		$ENV{'CPC_HOME'} = $deps->{'cpc_home'} if (!exists $ENV{'CPC_HOME'});

		if (!-d $ENV{'CPC_HOME'} || $ENV{'CPC_HOME'} eq '') {
		    $io->error("Environment variable CPC_HOME must be set pointing to CPC directory before running CPC.\n" .
			       'See CPC installation instructions (http://cpc.cbi.pku.edu.cn/docs/install_guide.jsp)');
		}
	    			    
		$io->execute_system_command("mkdir -p $cat_dir",
					    'Making output directory for ' . $io->file_basename($deps->{'cpc'}, 'suffix') . ' [ ' . $cat_dir . ' ]')
		    if (!-d $cat_dir);
				
		$io->execute_system_command("mkdir -p $cpc_work_dir") if (!-d $cpc_work_dir);

		my $cpc_cpu = 1;
		$cpc_cpu = $num_cpu_by_2 if (defined $num_cpu);
		
		if (!defined $skip_cpc_core) {
		    $io->execute_system_command('NUM_CPU=' . $cpc_cpu . ' ' . $deps->{'cpc'} . 
						' "' . $lncRNApipe_unique_seq_file . ' ' . $CPC_out_file . ' ' . $cpc_work_dir . ' ' . $cpc_work_dir . '" ' .
						$blastall_path->($deps->{'blastall'}),
						"Command call:\n" .
						"-------------\n" . 'NUM_CPU=' . $cpc_cpu . ' ' . $deps->{'cpc'} .
						' "' . $lncRNApipe_unique_seq_file . ' ' . $CPC_out_file . ' ' . $cpc_work_dir . ' ' . $cpc_work_dir . '" ' .
						$blastall_path->($deps->{'blastall'}));
		}
	    } if (defined $start_cpc);
	    
	    # To resume the remaining 2 modules, we need to check for CPC prediction file and adjust the file.

	    if (!defined $skip_cpc) {
		if (-e $CPC_out_file && -s $CPC_out_file) {
		    $io->c_time('Getting noncoding "only" transcripts [ from ' .
				$io->file_basename($unique_gtf, 'suffix') .
				" ] ...");
		    
		    $io->execute_system_command("grep noncoding $CPC_out_file | cut -f 1 | " .
						'while read trid; do grep -P "${trid//\./\\\.}\W" ' .
						"$unique_gtf; done &> $lncRNApipe_final_trs");
		    
		    $io->c_time('Creating FASTA for noncoding "only" transcripts [ from ' .
				$io->file_basename($unique_gtf, 'suffix') .
				" ] ...");
		    
		    $io->execute_system_command("grep noncoding $CPC_out_file | cut -f 1 | ".
						'while read trid; do grep -A 1 -P "${trid//\./\\\.}\W" ' .
						"$lncRNApipe_unique_seq_file; done &> $lncRNApipe_unique_seq_noncoding_file");
		    $lncRNApipe_unique_seq_file = $lncRNApipe_unique_seq_noncoding_file;
		}
		else {
		    $lncRNApipe_succ = 0;
		    $io->error("CPC prediction step seems to have failed... \n" .
			       "See --skip-cpc if you do not want to run CPC altogether\n\nBailing out!");
		}
	    } 
	    elsif (defined $skip_cpc) {
		$io->c_time('Creating a place holder [ --skip-cpc requested ] ...');
		$CPC_out_file = $lncRNApipe_final . '/.' . $io->file_basename($lncRNApipe_unique_seq_file) . '.ph';
		$io->execute_system_command('grep -oP ">.+?\\s+" ' . $lncRNApipe_unique_seq_file . " | sed -e 's/>//g' | sed -e 's/\\s+\$//g' " .
					    '| while read trid; do echo -e "${trid}\tnoncoding"; done > ' . $CPC_out_file);
		$io->execute_system_command("cp $unique_gtf $lncRNApipe_final_trs");
	    }

	    # Do RNAfold and draw entropy plot ############################################################

	    do {

		module_header($io, 'Module 6: Running RNAfold predictions [ on ' . $io->file_basename($lncRNApipe_unique_seq_file, 'suffix') . ' ] ...');
	
		$io->execute_system_command("mkdir -p $RNAfold_dir") if (!-d $RNAfold_dir);

		chdir $RNAfold_dir;

		$io->c_time('Making directory for RNAfold predictions [ ' . $io->file_basename($RNAfold_dir_final, 'suffix') . ' ] ...');
		$io->execute_system_command("mkdir -p $RNAfold_dir_final")			    
		    if (!-d $RNAfold_dir_final);

		$io->warning('CPC was not run! Will run RNAfold on all transcript models.') if (defined $skip_cpc);
	
		if (!defined $skip_rnafold_core) {
	
		    $io->c_time('Running RNAfold with -p flag [ on ' . $io->file_basename($lncRNApipe_unique_seq_file, 'suffix') .
				' ]. This may take very long time ...');
		    $io->execute_system_command($deps->{'RNAfold'} . 
						' -p ' . $start_rnafold . ' < ' . $lncRNApipe_unique_seq_file . ' &> ' . $RNAfold_mfe,
						"Command call:\n" .
						"-------------\n" . $deps->{'RNAfold'} .
						' -p ' . $start_rnafold . ' < ' . $lncRNApipe_unique_seq_file . ' &> ' . $RNAfold_mfe);
		}

		$io->c_time('Now generating RNAfold color plots ...');
		
		opendir (RNAfold_pl, $RNAfold_dir) || 
		    $io->error("Cannot open directory $RNAfold_dir to generate plots");
		
		push my @RNA_relPlots, grep {/\_dp\.ps$/} readdir RNAfold_pl;
				
		$io->error('Cannot find _dp PostScript files. Bailing out!')
		    if (scalar(@RNA_relPlots) == 0);

		foreach my $dp_ (@RNA_relPlots) {
		    $dp_ = $RNAfold_dir . '/' . $dp_;
		    my $ss_ = $RNAfold_dir . '/' . $io->file_basename($dp_, 'suffix');
		    $ss_ =~ s/\_dp\.ps$/\_ss\.ps/i;
		    my $ps = $ss_;
		    $ps =~ s/\_ss\.ps$/\_fstruct\.ps/;
		    my $pre_ps = $ps;
		    
		    $io->error('Cannot find corresponding _dp or _ss PostScript file(s) [ for transcript ' .
			       $io->file_basename($ps) . ' ] or any content in it. Bailing out!')
			if (!-e $ss_ || !-e $dp_ || !-s $ss_ || !-s $dp_);
				    
		    $ss_ = esc_tr_id($ss_);
		    $dp_ = esc_tr_id($dp_);
		    $ps = esc_tr_id($pre_ps);

		    $io->execute_system_command($deps->{'relplot'} . ' ' . $ss_ . ' ' . $dp_ . ' > ' . $ps,
						"Command call:\n" .
						"-------------\n" .
						$deps->{'relplot'} . ' ' . $ss_ . ' ' . $dp_ . ' > ' . $ps);
		    
		    $io->c_time('Cleaning up *_ss.ps and *_dp.ps files ...'), 
		    $io->execute_system_command("rm $dp_"),
		    $io->execute_system_command("rm $ss_")
			if ($rm_int_plots);
		    
		    $io->error('RNAfold prediction step failed. Cannot generate final PostScript file from _ss.ps and _dp.ps.' . "\nBailing out!") 
			if (!-e $pre_ps || !-s $pre_ps);

		    $io->execute_system_command("mv $ps $RNAfold_dir_final/.");

		}
		close RNAfold_pl;
		
		if (-e $CPC_out_file && -s $CPC_out_file) {
		    $io->c_time('Filtering transcripts [ ' . $io->file_basename($unique_gtf, 'suffix') . ' ] ...');    
		
		    $io->execute_system_command("grep noncoding $CPC_out_file | cut -f 1 | " .
						"while read trid; do find $RNAfold_dir -type f -name \${trid//\\./\\\\.}_fstruct.ps" .
						" -exec mv -t $RNAfold_dir_final" . "/\ {} \\;; done;");
		}

		$io->c_time('Moving around files and cleaning up directories ...');
		$io->execute_system_command("mv $RNAfold_dir/*.mfe $RNAfold_dir_final/.")
		    if (-e $RNAfold_mfe);

		# All possible checks complete, RNAfold has finished.
		$lncRNApipe_succ = 1;

	    } if (defined $start_rnafold);

	    # Start Infernal and edit final GTF  #############################################################

	    do {
		my $tbl = $final_infernal . '/' . $io->file_basename($lncRNApipe_unique_seq_file) . '.txt';
		my $raw = $final_infernal . '/' . $io->file_basename($lncRNApipe_unique_seq_file) . '.cmscan';

		$io->execute_system_command("mkdir -p $final_infernal")
		    if (!-d $final_infernal);		

		module_header($io, 'Module 7: Running cmscan from Infernal [ on ' . $io->file_basename($lncRNApipe_unique_seq_file, 'suffix') . ' ] ...');
	
		$deps->{'cmscan'} .= " --cpu $num_cpu_by_2" if (defined $num_cpu);

		if (!defined $skip_cmscan_core) {
		    
		    $io->execute_system_command($deps->{'cmscan'} . ' ' . $start_infernal . ' --tblout ' .
						$tbl . ' ' . $deps->{'rfam_cm'} . ' ' . $lncRNApipe_unique_seq_file . ' &> ' . $raw,
						"Command call:\n" .
						"-------------\n" . $deps->{'cmscan'} . ' ' . $start_infernal . ' --tblout ' .
						$tbl . ' ' . $deps->{'rfam_cm'} . ' ' . $lncRNApipe_unique_seq_file . ' &> ' . $raw);
		}

		$io->warning('CPC was not run! Will run Infernal on all transcript models.') if (defined $skip_cpc);
    
		if (-e $CPC_out_file && -s $CPC_out_file) {
		    $io->exist_sys_cmd(['sort', 'head', 'awk', 'uniq', 'mv', 'bc']);
		    $io->c_time('Updating final GTF file with Infernal annotations...');
		    
		    my $tmp_annot = $io->open_file('>', $lncRNApipe_final_trs . '.tmp');
		    unlink $lncRNApipe_allInfernal_trs if (-e $lncRNApipe_allInfernal_trs);
		    
		    my $infernal_run = $io->execute_get_sys_cmd_output("INF_GTF=$lncRNApipe_allInfernal_trs FINAL_GTF=$lncRNApipe_final_trs" .
								       " CM_TXT_OUT=$tbl CPC_TXT_OUT=$CPC_out_file COV=$inf_cov" .
								       ' ' . $deps->{'ann-inf'},
								       "Command call:\n" .
								       "-------------\n" .
								       "INF_GTF=$lncRNApipe_allInfernal_trs FINAL_GTF=$lncRNApipe_final_trs" .
								       " CM_TXT_OUT=$tbl CPC_TXT_OUT=$CPC_out_file COV=$inf_cov" .
                                                                       ' ' . $deps->{'ann-inf'});

		    $lncRNApipe_succ = pipeline_status($io, 'ERROR!', 'Cannot get annotation from Infernal run') if (!$infernal_run ||
														     $infernal_run =~ m/could not capture/i ||
														     !-e $tbl || !-e $raw || !-s $tbl || !-s $raw);
		    print $tmp_annot $infernal_run;
		    close $tmp_annot;
		    $io->c_time('Moving around files...');
		    $io->execute_system_command("mv $lncRNApipe_final_trs" . '.tmp ' . $lncRNApipe_final_trs);
		}
	      		    
	    } if (defined $start_infernal);

	    chomp (my $final_nc_tr_count = $io->execute_get_sys_cmd_output('grep -oP \'transcript_id\\s+\\".+?\\"\' ' . $lncRNApipe_final_trs .
									   ' | grep -oP \'\\".+\\"\' | sed -e \'s/\\"//g\' | sort -n | uniq | wc -l'));

	    if (!$final_nc_tr_count) {
		$io->execute_system_command( "cat $cat_dir/*.CPC.predict.txt | head", 
					     "Truncated output from  $cat_dir/*.CPC.predict.txt"); 
					    
		$io->error("Could not get final putative lncRNA count.\n" .
			   "Try re-running the final module. Take a look at --skip-blastall-core, --skip-rnafold-core and --skip-cmscan-core options." .
			   "\nThis may also mean that the putative lncRNAs that were used with CPC may have all been flagged as \"coding\".");
	    }
	    else {
		# Since we will have some transcripts by now that are filtered and passed all stages of the pipeline
		$lncRNApipe_succ = 1;
	    }

	    disp_ncRNA_counts($io, $lncRNApipe_final_trs, 'Final putative novel lncRNA count');
	    $cpu->finish(0, {'lncRNApipe_test' => $lncRNApipe_succ}) if (defined $num_cpu);

	    unlink $CPC_out_file if (defined $skip_cpc);
	}
	
	$cpu->wait_all_children if (defined $cpu);
	close lncRNApipe_fetch_seq_mod;
	
    } if (defined $start_cpc || defined $start_rnafold || defined $start_infernal || defined $skip_cpc);

    close $dep_tools_fh;
    close $pl_dep_paths_fh;

    # \x{2633}
    if ($lncRNApipe_succ) {
	$io->c_time("\x{2632}\x{2634} lncRNApipe Pipeline finished!");
    }
    else {
	$io->c_time("\x{2632}\x{2634} lncRNApipe Pipeline aborted(?)");
    }
    $io->end_timer($s_time);
}

################################################################################################
#                                                                                              #
#                                     Install the pipeline                                     #
#                                     ====================                                     #
#                                                                                              #
################################################################################################

sub setup_lncRNApipe {
    chomp(my $install_dir = shift);
    chomp(my $dl_util = `which wget 2>&1`);

    binmode(STDOUT, ":utf8");
    print "\n\x{1f473}  Hi! I am Vayu and I will attempt to install lncRNApipe for you.";
    print "\nChecking for platform independent prerequisites on UNIX based machines...\n";

    # Wget or Curl
    if ($dl_util !~ m/.+?wget$/i) {
	chomp ($dl_util = `which curl 2>&1`);
	succ_or_fail(0, 'Curl|Wget')
	    if ($dl_util !~ m/.+?curl$/i);
	succ_or_fail(1, 'curl');
        $dl_util .= ' -LkO ';
    }
    elsif ($dl_util =~ m/.+?(wget)$/i) {
	succ_or_fail(1, $1);
	$dl_util .= ' --no-check-certificate ';
    }
    else{
	succ_or_fail(0, 'Curl|Wget');
    }

    # Unzip
    check_util('unzip');

    # Make
    check_util('make');

    # Echo
    check_util('echo');

    # Uname
    check_util('uname');

    # bunzip2
    check_util('bunzip2');

    # gunzip
    check_util('gunzip');

    # rm
    check_util('rm');

    # cp
    check_util('cp');

    # mkdir
    check_util('mkdir');

    # touch
    check_util('touch');
    
    # mv
    check_util('mv');

    # find
    check_util('find');

    # bc
    check_util('bc');

    print "\n\nWe definetely need the custom module, IO::Routine.\n" .
      "We also need good grace of cpanm to install remaining modules.\n\n" .
	"Attempting to fetch IO::Routine from Perl-for-Bioinformatics repository on github...\n";

    if (-d "$install_dir/.build") {
	system("rm -rf $install_dir/.build");
    }

    if (-d '.build') {
	system("rm -rf .build");
    }

    system("mkdir $install_dir/.build");
    succ_or_fail(0, '.build dir') if (!-d "$install_dir/.build");

    system("rm $install_dir/master.zip") if (-e "$install_dir/master.zip");

    my $dl_util_slave = '';
    if ($dl_util =~ m/.*?wget\s+/i) {
      $dl_util_slave = $dl_util . '-O master.zip';
    }
    else {
      $dl_util_slave = $dl_util;
    }
    
    print "\n$dl_util_slave https://github.com/biocoder/Perl-for-Bioinformatics/archive/master.zip\n\n";
    system("$dl_util_slave https://github.com/biocoder/Perl-for-Bioinformatics/archive/master.zip > /dev/null 2>&1");
    system("mv $install_dir/master.zip $install_dir/.build/master.zip");

    succ_or_fail(0, 'master.zip')
	if (!-e "$install_dir/.build/master.zip");

    print "Inflating master.zip...\n";
    print "\nunzip -d $install_dir/.build $install_dir/.build/master.zip\n";
    system("unzip -d $install_dir/.build $install_dir/.build/master.zip  > /dev/null 2>&1");

    # Use cpanm to install Set::IntervalTree
    system("$dl_util https://raw.githubusercontent.com/miyagawa/cpanminus/master/cpanm  > /dev/null 2>&1");
    succ_or_fail(0, 'cpanm') if (!-e 'cpanm');
    system("mv cpanm $install_dir/.build");
    system("chmod 755 $install_dir/.build/cpanm");

    my $pm_install_dir = $install_dir;
    if (defined $local_lib) {
	chop $local_lib if ($local_lib =~ m/\/$/);
	$pm_install_dir = $local_lib;
    }

    cpanm_status($pm_install_dir, 'Module::Load', $install_dir);
    cpanm_status($pm_install_dir, 'Set::IntervalTree', $install_dir);
    cpanm_status($pm_install_dir, 'LWP::Simple', $install_dir);
    cpanm_status($pm_install_dir, 'XML::XPath', $install_dir);
    cpanm_status($pm_install_dir, 'XML::XPath::XMLParser', $install_dir);
    cpanm_status($pm_install_dir, 'Parallel::ForkManager', $install_dir);
    cpanm_status($pm_install_dir, 'MIME::Lite', $install_dir);
    cpanm_status($pm_install_dir, 'Email::Valid', $install_dir);

    if (!ask_user("\n\nDo you already have BioPerl installed?")) {
	print "\nOK ... This may take a while ... Grab a snack or coffee!\n";
	cpanm_status($pm_install_dir, 'BioPerl', $install_dir);
    }


    print "\nThank you cpanm!\n";

    my $pm_path = "$pm_install_dir/PERLLIBS:$pm_install_dir/PERLLIBS/lib:$pm_install_dir/PERLLIBS/lib/perl5";

    print "\nInstalling IO::Routine...\n\n";
    system("mkdir $pm_install_dir/PERLLIBS") if (!-d "$pm_install_dir/PERLLIBS");
    my $custom_pm_log = `cd $install_dir/.build/Perl-for-Bioinformatics-master/IO-Routine;perl Makefile.PL PREFIX=$pm_install_dir/PERLLIBS LIB=$pm_install_dir/PERLLIBS/lib && make && make test && make install; 2>&1`;

    if ($custom_pm_log =~ m/fail|error|cannot/i) {
	succ_or_fail(0, 'IO::Routine');
    }
    else {
	succ_or_fail(1, 'IO::Routine');
    }

    require "$pm_install_dir/PERLLIBS/lib/IO/Routine.pm";
    my $io = IO::Routine->new($help, $quiet);

    chomp(my $lncRNApipe_root = dirname(abs_path($0)));

    if (!-d "$lncRNApipe_root/.lncRNApipe.depbin") {
	print "\n\nMoving required tools and scripts to current install path...\n";
	system("mkdir $install_dir/.lncRNApipe.depbin");
	system("mv $install_dir/.build/Perl-for-Bioinformatics-master/NGS-Utils/.lncRNApipe.depbin/* $install_dir/.lncRNApipe.depbin/.");
	system("mv $install_dir/.build/Perl-for-Bioinformatics-master/NGS-Utils/* $install_dir/.");
    }
    
    print "\n\nCleaning up build directory...\n";
    system("rm -rf $install_dir/.build");

    $io->execute_system_command(0, 'Detecting system architecture...');

    my $sys_arch_info = $io->execute_get_sys_cmd_output('uname -a');

    $io->error('This is not a UNIX based machine ... Aborting installation!')
	if (!$sys_arch_info ||
	    $sys_arch_info !~ m/linux|darwin/i);

    $io->error('It is not a 64-bit machine!' . 
	       "\nIt is your responsibility to make sure that you have the following tools" .
	       " installed for your system architecture and also must be found in \$PATH:\n" .
	       "\nblastall, gtfToGenePred, cuffcompare and RNAfold\n\nBailing out!") 
	if ($sys_arch_info !~ m/x86\_64/i);


    my $sys_arch = '';
    if ($sys_arch_info =~ m/darwin/i) {
	$sys_arch = 'darwin';
	$io->execute_system_command(0, 'Skipping version requirement check for system level commands [ tar ], [ cut ] and [ wc ].' .
				    "\nFreeBSD's tools does not provide version numbers (?)");
	check_util('tar');
	check_util('cut');
	check_util('wc');

	# Darwin special case...
	chomp(my $darwin_sed = `sed --version 2>&1`);
	chomp(my $darwin_grep = `grep --version 2>&1`);

	if ($darwin_sed !~ m/gnu/i || $darwin_grep !~ m/gnu/i) {
	    print "\n";
	    $io->warning("We need GNU's grep and sed instead of FreeBSD's.");
	    chomp(my $homebrew = `which brew 2>&1`);
	    print "Looking to see if you have homebrew to install GNU tools [ grep and sed ] ...\n";
	    
	    if ($homebrew =~ m/.*brew$/i) {
		succ_or_fail(1, 'homebrew');
	    }
	    else {
		$io->error('Please install homebrew. See installation instructions at http://brew.sh/' . 
			   "\nThen, do the following:\n\nbrew update\nbrew tap homebrew/dupes\n" .
			   "brew install grep\nbrew install gnu-sed\n" .
			   "ln -s /usr/local/bin/ggrep /usr/local/bin/grep\n" .
			   "ln -s /usr/local/bin/gsed /usr/local/bin/sed\n\n" .
			   "Then, rerun the setup procedure to successfully install lncRNApipe pipeline.");
	    }
	    print "\n\nPlease execute the following commands in order and" .
		" rerun the setup procedure to successfully install lncRNApipe pipeline.\n";
	    $io->execute_system_command(0,
					"brew update\nbrew tap homebrew/dupes\n" .
					"brew install grep\nbrew install gnu-sed\n" .
					"ln -s /usr/local/bin/ggrep /usr/local/bin/grep\n" .
					"ln -s /usr/local/bin/gsed /usr/local/bin/sed");
	    print "Skipping version check for GNU grep and GNU sed...\n";
	    succ_or_fail(2, 'sed');
	    succ_or_fail(2, 'grep');
	    print "\n\nAborting setup...\n\n";
	    exit;
	}
	else {
	    $io->check_sys_level_cmds(['grep', 'sed'],
				      ['2.6.3', '4.2.1']);
	    succ_or_fail(1, 'grep');
	    succ_or_fail(1, 'sed');
	}
    }
    else {
	$sys_arch = 'linux';
	 $io->check_sys_level_cmds(['grep', 'sed', 'tar', 'cut', 'wc'],
				   ['2.6.3', '4.2.1', '0', '8', '8']);

	succ_or_fail(1, 'grep');
	succ_or_fail(1, 'sed');
	succ_or_fail(1, 'tar');
	succ_or_fail(1, 'cut');
	succ_or_fail(1, 'wc');
    }

    print "\n\nSetting up PERL5LIB paths...\n";
    unlink "$USER_HOME/.lncRNApipe.PERLLIBS" if (-e "$USER_HOME/.lncRNApipe.PERLLIBS");
    my $pl_dep_fh = $io->open_file('>', "$USER_HOME/.lncRNApipe.PERLLIBS");
    print $pl_dep_fh $pm_path;

    $io->execute_system_command(0,
				"\nChecking for lncRNApipe pipeline tool dependencies..." .
				"\n\nWriting tool dependency chain to $USER_HOME/.lncRNApipe.depconf");

    my $dep_fh = $io->open_file('>', "$USER_HOME/.lncRNApipe.depconf");

    # We will figure out if user is installing at a different location other than
    # from the cloned repo.
    $install_dir = $lncRNApipe_root if (-d "$lncRNApipe_root/.lncRNApipe.depbin" && $lncRNApipe_root !~ m/^\./);

    if (-e "$install_dir/.lncRNApipe.depbin/Rfam.cm.1_1.bz2" &&
	-e "$install_dir/.lncRNApipe.depbin/Rfam.cm.1_1.i1m.bz2") {
	print "Decompressing Rfam CM files...\n\n";
	system("bunzip2 $install_dir/.lncRNApipe.depbin/Rfam.cm.1_1.bz2");
	system("bunzip2 $install_dir/.lncRNApipe.depbin/Rfam.cm.1_1.i1m.bz2");
    }

    print $dep_fh check_bio_util('cuffcompare', $install_dir, $sys_arch), "\n";
    print $dep_fh check_bio_util('blastall', $install_dir, $sys_arch), "\n";
    print $dep_fh check_bio_util('RNAfold', $install_dir, $sys_arch), "\n";
    print $dep_fh check_bio_util('gtfToGenePred', $install_dir, $sys_arch), "\n";
    print $dep_fh check_bio_util('cmscan', $install_dir, $sys_arch), "\n";
    print $dep_fh check_bio_util('formatdb', $install_dir, $sys_arch), "\n";
    print $dep_fh "$install_dir/" . $CAT_SC, "\n" if (check_native("$install_dir/" . $CAT_SC));
    print $dep_fh "$install_dir/" . $FETCH_SC, "\n" if (check_native("$install_dir/" . $FETCH_SC));
    print $dep_fh "$install_dir/" . $GET_UQ_SC, "\n" if (check_native("$install_dir/" . $GET_UQ_SC));
    print $dep_fh "$install_dir/" . $RELPLOT_MOD_SC, "\n" if (check_native("$install_dir/" . $RELPLOT_MOD_SC));
    print $dep_fh "$install_dir/" . $CPC_SH_SC, "\n" if (check_native("$install_dir/" . $CPC_SH_SC));
    print $dep_fh "$install_dir/" . $ANN_INF_SC, "\n" if (check_native("$install_dir/" . $ANN_INF_SC));
    print $dep_fh "$install_dir/" . $SEND_MAIL_SC, "\n" if (check_native("$install_dir/" . $SEND_MAIL_SC));

    if (-e "$install_dir/.lncRNApipe.depbin/Rfam.cm.1_1" &&
	-s "$install_dir/.lncRNApipe.depbin/Rfam.cm.1_1") {
	print $dep_fh "$install_dir/.lncRNApipe.depbin/Rfam.cm.1_1\n";
	succ_or_fail(1, "Rfam.cm.1_1");
    }
    else {
	succ_or_fail(0, "Rfam.cm.1_1");
    }

    if (!exists $ENV{'CPC_HOME'}) {
	if (ask_user("\n\nDo you want me to attempt to install CPC?")) {
	    	setup_cpc($dl_util_slave, $dep_fh);
		succ_or_fail(1, 'CPC');
	}
	else {
	    succ_or_fail(2, 'CPC');
	}
	#print "\n\nSeems like Coding Potential Calculator (CPC) is not installed. This may not be a problem if you do not intend to use CPC.";
	#print "\nSee CPC installation instructions (http://cpc.cbi.pku.edu.cn/docs/install_guide.jsp) if you would like to run CPC module with lncRNApipe.";
    }
    close $dep_fh;

    $io->execute_system_command(0,
				"\n\nSetup complete. See \"perl lncRNApipe --h\" to run the pipeline with options.");
    # For build check
    exit;
}


# Try and install CPC into current install path. 
sub setup_cpc {
    my $dl_util = shift;
    my $fh = shift;

    return if (exists $ENV{'CPC_HOME'});

    my $depconf_fh = read_conf('/.lncRNApipe.depconf');
 
    my $curr_dir = cwd . '/.cpc';
    system("rm -rf $curr_dir") if (-d $curr_dir);
    
    print $fh "cpcHome:$curr_dir\n";
    $ENV{'CPC_HOME'} = $curr_dir;
    
    get_deps($depconf_fh);

    print "\n\nAttempting to setup CPC [ github.com/biocoder/cpc ]. \nThis may take a while. Stay put...\n\n";
    system("$dl_util https://github.com/biocoder/cpc/archive/master.zip > /dev/null 2>&1");
    system("unzip -d $curr_dir master.zip > /dev/null 2>&1");
    system("rm master.zip");
    system("mv $curr_dir/cpc-master/* $curr_dir/.");
    rmdir "$curr_dir/cpc-master";

    chdir "$curr_dir/cpc-master/libs/libsvm/libsvm-2.81";
    system("make clean && make > /dev/null 2>&1");

    chdir "$curr_dir/cpc-master/libs/estate";
    system("make clean && make > /dev/null 2>&1");

    my $format_db = 'uniref90.fasta.gz';
    if (!ask_user("\nIf you already have uniref90.fasta file somewhere, then there is no point in wasting time to download it again.\n" .
		  "Do you already have uniref90.fasta somewhere on this system?")) {
	
	print "\nAttempting protein database download...\n\n";
	$dl_util =~ s/master\.zip/uniref90\.fasta\.gz/;
	system("$dl_util ftp://ftp.ebi.ac.uk/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz");
    }
    else {
	$format_db = ask_user("\nEnter full path to uniref90.fasta: ");
	succ_or_fail(0, "CPC [ ERROR: $! ]")
	    if (!-r $format_db || !-e $format_db || !-s $format_db);
    }

    print "\nDecompressing uniref90.fasta.gz...\n",
    system("gunzip uniref90.fasta.gz") if ($format_db =~ m/\.gz$/);
	
    if (!-e "$curr_dir/data/formatdb.log" ||
	!-s "$curr_dir/data/formatdb.log" ||
	!-e "$curr_dir/data/prot_db.pal" ||
	!-s "$curr_dir/data/prot_db.pal") {
	print "\nFormatting database...\n";
	chdir "$curr_dir/data";
	system($deps->{'formatdb'} . " -i $format_db -p T -n prot_db > /dev/null 2>&1");
    }
    
    return;
}


# Check for pipeline scripts.
sub check_native {
    my $nat_sc = shift;
    
    (my $sc) = ($nat_sc =~ m/.*\/(.+?\.pl|.+?\.sh)$/);

    if (-e $nat_sc) {
	succ_or_fail(1, $sc);
	return 1;
    }
    else {
	succ_or_fail(0, $sc);
    }
    return 0;
}

# Check for dependencies
sub check_bio_util {
    my $cmd = shift;
    my $install_dir = shift;
    my $sys_arch = shift;

    chomp (my $biocmd4arch = `which $cmd 2>&1`);
    my $sandboxed4arch = $sys_arch . '_' . $cmd;

    if ($biocmd4arch !~ m/.*$cmd$/i &&
        -d "$install_dir/.lncRNApipe.depbin" &&
        -e "$install_dir/.lncRNApipe.depbin/$sandboxed4arch") {
        succ_or_fail(1, $cmd);
	$biocmd4arch = "$install_dir/.lncRNApipe.depbin/$sandboxed4arch";
    }
    elsif ($biocmd4arch =~ m/.*$cmd$/i) {
	succ_or_fail(1, $cmd);
    }
    else {
	succ_or_fail(0, $cmd);
    }
    return $biocmd4arch;
}

# Check for system utilities
sub check_util {
    my $cmd = shift;
    chomp (my $util = `which $cmd 2>&1`);

    if ($util !~ m/.+?$cmd/i) {
	succ_or_fail(0, $cmd);
	print "\nAborting...\n\n";
	exit;
    }
    else {
        succ_or_fail(1, $cmd);
    }
    return $util;
}

# Print success or failure of lncRNApipe pipeline modules
sub pipeline_status {
    my $self = shift;
    my $status = shift;
    my $msg = shift;

    if ($status =~ m/ERROR|STDERR/i) {
	if ($msg) {
	    print "\nERROR!\n------\n$msg\n";
	}
	else {
	    print "\n$status\n";
	}
	$self->c_time("\x{2632}\x{2634} lncRNApipe Pipeline aborted(?)\n\n");
	exit;
    }
    elsif ($status !~ m/STDERR/i) {
	print $status;
    }
    return 1;
}

# Print success or failure with character length adjustments
sub succ_or_fail {
    my $code = shift;
    my $msg = shift;

    my $char_white_space = 70 - length($msg);
    $char_white_space = 5 if ($char_white_space < 0);

    print "\n";
    printf("%s%*s%s", $msg, $char_white_space, '... ', BOLD GREEN . 'OK' . RESET)
	if ($code == 1);

    if ($code == 0) {
	printf("%s%*s%s", $msg, $char_white_space, '... ',  BOLD RED . 'FAIL' . RESET);
	print "\n\n";
	exit;
    }

    printf("%s%*s%s", $msg, $char_white_space, '... ', BOLD YELLOW . 'SKIP' . RESET)
         if ($code == 2);

    printf("%s%*s%s", $msg, $char_white_space, '... ', BOLD YELLOW . 'WARN' . RESET)
	if ($code == 3);

    return;
}

# Format module header to easily spot the module start in log output
sub module_header {
    my $c_io = shift;
    my $heading = shift;
    
    my $char_hash_len = 73 - length($heading);
    $char_hash_len = 0 if ($char_hash_len < 0);
    my $char_hash_fixed = '#' x 27;
    my $char_hash_var = '#' x $char_hash_len;

    my $heading_string = sprintf("%s%s%s", $char_hash_fixed, ' ' . $heading . ' ', $char_hash_var);
    
    $c_io->c_time($heading_string);
    return;
}

# Source the environment
sub source_env {
    chomp(my $pl_dep_paths_fh = shift);
    chomp(my $pl_paths = do { local $/; <$pl_dep_paths_fh> });

    my $perl_inc_string = ' ';
    
    foreach my $pl_path (split(/\:/, $pl_paths)) {
	$perl_inc_string .= '-I' . $pl_path . ' ';
	push @INC, $pl_path;
    }
    return $perl_inc_string;
}

# Report cpanm status
sub cpanm_status {
    my $install_dir = shift;
    my $pm = shift;
    my $cpanm_bin = shift;

    my $cpanm_log = `$cpanm_bin/.build/cpanm -f -l $install_dir/PERLLIBS $pm 2>&1`;
    
    if ($cpanm_log !~ m/.*?success.*/i && $pm =~ m/xml/i) {
	print STDOUT "\n\nYou need to install XML parser C libraries.\n\n\t\* On Ubuntu / Debian based Linux distributions, as root user, do:\n\n\t\tapt-get install libexpat1 libexpat1-dev\n\n\t\* On RedHat / Fedora / CentOS based Linux distributions, as root user do:\n\n\t\tyum install expat expat-devel\n\n";
	succ_or_fail(0, $pm);
    }

    succ_or_fail(0, $pm) if ($cpanm_log !~ m/.*?success.*/i);
    succ_or_fail(1, $pm);
    return;
}

# Report transcript files and size
sub get_num_tr_files {
    my $cmd_args = shift;
    my @args = split/\s+/, $cmd_args;
   
    my $last_opt_elem = my $tr_size = 0;
    my $tr_files = '';
    
    for (0 .. $#args) {
	$last_opt_elem = $_ + 2 if ($args[$_] =~ m/^\-/);
    }
    $tr_size = $#args - $last_opt_elem;
    
    $tr_files = join(' ', @args[$last_opt_elem .. $#args]);
    return($tr_files, $tr_size);
}

# Get number of transcripts
sub get_num_trs {
    my $self = shift;
    my $file = shift;
    $self->exist_sys_cmd(['wc', 'sort', 'uniq']);
    return $self->execute_get_sys_cmd_output('grep -oP \'transcript_id \\".+?\\"\' ' . $file .
					     ' | awk \'{print $2}\' | sed -e \'s%\"%%g\' | sort -n | uniq | wc -l');
}

# Get cuffcompare's annotation file
sub get_annot {
    my $cmd_args = shift;
    my $annot_file = '';
    
    ($annot_file) = ($cmd_args =~ m/.*?-r\s+(.+?)\s+.*/i);
    
    return $annot_file;
}

# Ask for user input
sub ask_user {
 my $msg = shift;
 my $ans;

 while (1) {
     print $msg, ' [ y/n ]: ';
     
     eval {
	 local $SIG{ALRM} = sub { die "DidNotAnswer" };
	 alarm 60;
	 chomp($ans = <STDIN>);
     };
     
     succ_or_fail(0, "\n\nAborted! [ No answer for 60 seconds. ] ??") if ($@ =~ m/DidNotAnswer/i);
     
     if ($ans !~ m/^y|n|fasta|fa$/i) {
	 print "\nInvalid answer!\n";
	 next;
     }
     else {
	 alarm 0;
	 last;
     }
 }
 
 return 1 if ($ans =~ m/^y$/i);
 return 0 if ($ans =~ m/^n$/i);
 return $ans if ($ans =~ m/fasta|fa$/);
}

# Get and verify dependencies
sub get_deps {
    my $dep_tools_fh = shift;
    my $dep_tools = do {local $/; <$dep_tools_fh>};
    ($deps->{'cuffcompare'}) = ($dep_tools =~ m/(.+?cuffcompare)/i);
    ($deps->{'blastall'}) = ($dep_tools =~ m/(.+?blastall)/i);
    ($deps->{'RNAfold'}) = ($dep_tools =~ m/(.+?rnafold)/i);
    ($deps->{'bin-gtfToGenePred'}) = ($dep_tools =~ m/(.+?gtfToGenePred)/i);
    ($deps->{'cmscan'}) = ($dep_tools =~ m/(.+?cmscan)/i);
    ($deps->{'formatdb'}) = ($dep_tools =~ m/(.+?formatdb)/i);
    ($deps->{'rfam_cm'}) = ($dep_tools =~ m/(.+?Rfam\.cm\.1\_1)/i);
    ($deps->{'cat'}) = ($dep_tools =~ m/(.+?$CAT_SC)/i);
    ($deps->{'fetch'}) = ($dep_tools =~ m/(.+?$FETCH_SC)/i);
    ($deps->{'get'}) = ($dep_tools =~ m/(.+?$GET_UQ_SC)/i);
    ($deps->{'cpc'}) = ($dep_tools =~ m/(.+?$CPC_SH_SC)/i);
    ($deps->{'ann-inf'}) = ($dep_tools =~ m/(.+?$ANN_INF_SC)/i);
    ($deps->{'relplot'}) = ($dep_tools =~ m/(.+?$RELPLOT_MOD_SC)/i);
    ($deps->{'email'}) = ($dep_tools =~ m/(.+?$SEND_MAIL_SC)/i);
    ($deps->{'cpc_home'}) = ($dep_tools =~ m/cpcHome\:(.+)/i) if (!exists $ENV{'CPC_HOME'}); 

    foreach my $dep (keys %$deps) {
	die "\n\nCannot find required dependency [ " . $dep . " ]. May be reinstalling lncRNApipe will resolve the issue?\n\n"
	    if (!$deps->{$dep} || !-e $deps->{$dep});
    }
    return;
}

# Read lncRNApipe configuration
sub read_conf {
    my $conf = shift;
    open(my $conf_fh, '<', $USER_HOME . $conf) ||
        die "\nCannot open $USER_HOME$conf" . ".\nMay be setup was unsuccessful (?): $!\!\n" .
        "\nRun \"perl lncRNApipe -setup\" to try to setup lncRNApipe pipeline.\n\n";
    return $conf_fh;
}

# Escape tr ids
sub esc_tr_id {
    my $tr_id = shift;

    $tr_id =~ s/\|/\\\|/;
    $tr_id =~ s/\:/\\\:/;

    return $tr_id;
}

# Load parallel fork module
sub do_parallel {
    my $num_cpu = shift;
    require Parallel::ForkManager;
    Parallel::ForkManager->import();
    my $cpu = Parallel::ForkManager->new($num_cpu);
    $cpu->set_max_procs($num_cpu);
    return $cpu;
}

# Display ncRNA category counts
sub disp_ncRNA_counts {
    my $c_io = shift;
    my $file = shift;
    my $total_desc = shift;

    chomp (my $num_lincs = $c_io->execute_get_sys_cmd_output("grep -iP '\ttranscript\t.+?lncRNA_type.+?lincrna.*?Infernal' $file | wc -l"));
    $num_lincs = 0 if ($num_lincs =~ m/^could not capture/i);

    chomp (my $num_concs = $c_io->execute_get_sys_cmd_output("grep -iP '\ttranscript\t.+?lncRNA_type.+?conc.*?.*?Infernal' $file | wc -l"));
    $num_concs = 0 if ($num_concs =~ m/^could not capture/i);
    
    chomp (my $num_poncs = $c_io->execute_get_sys_cmd_output("grep -iP '\ttranscript\t.+?lncRNA_type.+?ponc.*?.*?Infernal' $file | wc -l"));
    $num_poncs = 0 if ($num_poncs =~ m/^could not capture/i);    

    chomp (my $num_incs = $c_io->execute_get_sys_cmd_output("grep -iP '\ttranscript\t.+?lncRNA_type.+?\"inc.+?intronic.*?.*?Infernal' $file | wc -l"));
    $num_incs = 0 if ($num_incs =~ m/^could not capture/i);

    chomp (my $num_ex_ov = $c_io->execute_get_sys_cmd_output("grep -iP '\ttranscript\t.+?lncRNA_type.+?exonic\\s+overlap.*?.*?Infernal' $file | wc -l"));
    $num_ex_ov = 0 if ($num_ex_ov =~ m/^could not capture/i);

    chomp (my $num_inf_ann = $c_io->execute_get_sys_cmd_output("grep -viP '\ttranscript\t.+?No Annotation' $file | grep -iP '\ttranscript\t' | wc -l"));
    $num_inf_ann = 0 if ($num_inf_ann =~ m/^could not capture/i);

    print("\n\nlncRNA Summary [ " . $c_io->file_basename($file, 'suffix') . " ] :\n" .
	  "---------------------------------------------------------------------------------------\n" .
	  "LincRNAs: $num_lincs\n" .
	  "Intronic overlaps - Concs: $num_concs\n" .
	  "Intronic overlaps - Poncs: $num_poncs\n" .
	  "Intronic overlaps - Incs: $num_incs\n" .
	  "Exonic overlaps: $num_ex_ov\n" .
	  $total_desc . ': ' . ($num_lincs +
				$num_concs +
				$num_incs +
				$num_ex_ov +
				$num_poncs) .
	  "\nInfernal matches: $num_inf_ann\n" .
	  "\n\n" );
    
    
    return;
    
}

# Check for new release on github
sub check_for_updates { 
    return if (defined $no_update_check);
    my $io = run_lncRNApipe('version');

    $io->warning('Checking for updates on github.com/biocoder ...', 'INFO!!');
    my $newest_release = $io->execute_get_sys_cmd_output("curl -Lks https://api.github.com/repos/biocoder/Perl-for-Bioinformatics/releases | grep 'tag_name'");
    my $newest_release_date = $io->execute_get_sys_cmd_output("curl -Lks https://api.github.com/repos/biocoder/Perl-for-Bioinformatics/releases | grep 'created_at' | head -n 1");
    
    if ($newest_release =~ m/^could not/i) {
	$io->warning('Unable to contact github.com for updates ... May be we have reached the X-rate-limit! :-(' .
		     "\nUse --d option if you wish to disable update check.");
	return;
    }
    
    $newest_release =~ s/\"|\,|\v//g;
    $newest_release_date =~ s/\"|\,|\v|Z//g;
    $newest_release_date =~ s/T/ /;
    ($newest_release) = ($newest_release =~ m/.+?name.+?(\d.*)/);
    ($newest_release_date) = ($newest_release_date =~ m/\:(.+)/);
    $newest_release_date = $io->strip_leading_and_trailing_spaces($newest_release_date);
       
    my $is_upToDate = $io->check_sys_level_cmds(['lncRNApipe'],
						["$newest_release"], 'warn'); 
    
    $io->warning('You are using the latest version of the software which is v' . $newest_release . " [ Release date: $newest_release_date  ]." .
		 "\nUse --d option if you wish to disable update check.", 'OK!') if (!$is_upToDate);

    chomp(my $last_commit_msg = $io->execute_get_sys_cmd_output("curl -Lks https://api.github.com/repos/biocoder/Perl-for-Bioinformatics/commits?per_page=1 | grep 'message'"));
    chomp(my $last_commit_date = $io->execute_get_sys_cmd_output("curl -Lks https://api.github.com/repos/biocoder/Perl-for-Bioinformatics/commits?per_page=1 | grep 'date' | head -n 1"));

    if ($last_commit_msg =~ m/^could not/i) {
	$io->warning('Unable to contact github.com for latest commit message ... May be we have reached the X-rate-limit! :-(' .
		     "\nUse --d option if you wish to disable update check.");
	return;
    }

    $last_commit_msg =~ s/^\s+|\"|\,$//g;
    $last_commit_msg =~ s/message:\s*/Last Commit: /;
    $last_commit_msg .= '.' if ($last_commit_msg !~ m/\.$/);
    $last_commit_date =~ s/\"|\,|\v|Z//g;
    $last_commit_date =~ s/T/ /g;
    ($last_commit_date) = ($last_commit_date =~ m/\:(.+)/);
    $last_commit_date = $io->strip_leading_and_trailing_spaces($last_commit_date);

    if ($VERSION !~ m/$newest_release/) {
	$io->warning($last_commit_msg .
		     "\nVersion: v$newest_release\nCommit date: $last_commit_date",
		     'NEW!');
    }    
    else {
	$io->warning($last_commit_msg, 'COMMIT');
    }
    return;
}

# Get lncRNApipe pipeline status code after child exits
sub get_pipeline_status_code {
    my $cpu = shift;
    return $cpu->run_on_finish( 
	sub {
	    my ($pid, $exit_code, $p_identifier, $exit_sig, $core_dump, $stat) = @_;
	    return $stat->{'lncRNApipe_succ'};
	});
}


__END__

=head1 NAME

 _            ____  _   _    _          _            
| |_ __   ___|  _ \| \ | |  / \   _ __ (_)_ __   ___ 
| | '_ \ / __| |_) |  \| | / _ \ | '_ \| | '_ \ / _ \
| | | | | (__|  _ <| |\  |/ ___ \| |_) | | |_) |  __/
|_|_| |_|\___|_| \_|_| \_/_/   \_| .__/|_| .__/ \___|
                                 |_|     |_|         

=head1 SYNOPSIS

lncRNApipe is a pipeline to extract putative novel lncRNAs ab initio, given a list of transcripts in GTF format assembled from deep sequencing data (ex: RNA-Seq).

Complete Description: 

 perldoc lncRNApipe

Examples:

 perl lncRNApipe --h

 perl lncRNApipe --run /data/lncRNApipe/ \
                 --cuffcmp '-r annotation.gtf -s genome.fa transcripts1.gtf transcripts2.gtf' \
                 --cat-ncRNAs '-sample-names "M1,M2" -ov 80 -fpkm 2 -len 200 -max-len 10000 -min-exons 1 -antisense" \
                 --get-uq-feat '-sf known_lncRNAs.bed' \
                 --fetch-seq '-db mm10' \
                 --cpc \
                 --rna \
                 --inf &> lncRNApipe.run.log

=head1 DESCRIPTION

This pipeline script will bind together the functionality of the tools / scripts: cuffcompare, categorize_ncRNAs.pl, get_unique_features.pl, fetch_seq_from_ucsc.pl, RNAfold, Infernal and Coding Potential Calculator (CPC). Transcriptome construction tools such as Cufflinks produces a set of assembled transcripts in GTF format. lncRNApipe uses this data in addition to known gene annotation to extract putative lncRNAs constructed by the ab initio assemblers. The pipeline relies on the FPKM / RPKM values generated by these assemblers to assess the confidence of the constructed de novo transcripts and validates it against the known reference gene and non coding RNA information to identify putative novel lncRNAs. In brief, the pipeline steps are as follows:

=over 5

=item 1. Cuffcompare (lncRNApipe Option: --cuff or --cuffcompare)

The transcript assembly can be compared to known annotation of choice to classify them into different class codes (http://cufflinks.cbcb.umd.edu/manual.html#class_codes) using Cuffcompare. The assembled transcripts should be in GTF format. Cufflinks generates the output files in GTF format. If you are using other software such as Scripture, you can convert the output file in BED format into GTF using Scripture as:

 java -jar /path/to/scripture.jar -task toGFF -cufflinks -in scripture.out.bed -source SCRIPTURE -out scripture.out.gtf

To get help documentation about this module, run:

 perl lncRNApipe --h cuff

=item 2. categorize_ncRNAs.pl (lncRNApipe Option: --cat or --cat-ncRNAs)
    
lncRNApipe uses the tracking file (*.tracking) produced by Cuffcompare, annotation data in gene prediction format and a list of supplied transcripts in GTF format to produce and categorize lncRNAs into different classes as mentioned in the paper: http://genome.cshlp.org/content/22/3/577.full. In brief, the lncRNAs are classified into 5 categories: Long intergenic lncRNAs (LincRNAs), Intronic contained lncRNAs (Incs), Partially overlapping lncRNAs (Poncs), Completely overlapping lncRNAs (Concs) and Exonic overlaps (LncRNAs with sense or antisense overlap with reference exon).

To get help documentation about this module, run:

 perl lncRNApipe --h cat

=item 3. get_unique_features.pl (lncRNApipe Option: --get or --get--uq-feat)
    
It then compares the putative list with supplied known lncRNAs in BED format to get features that do not overlap any known lncRNAs. To extract known lncRNAs from your assembled transcripts for downstream analysis, include --known with --get option of lncRNAscan (--known option with lncRNApipe is experimental).

To get help documentation about this module, run:

 perl lncRNApipe --h get

=item 4. fetch_seq_from_ucsc.pl (lncRNApipe Option: --fetch or --fetch-seq)
    
This list is then used to fetch DNA sequence of those transcript sequences to determine their coding potential using CPC.
    
To get help documentation about this module, run:
    
 perl lncRNApipe --h fetch

=item 5. CPC.sh (lncRNApipe Option --cpc or --cpc)
    
In this step, the fetched FASTA sequences are used to determine the coding potential and those that are flagged as "noncoding" are used to create the final list of high confidence lncRNAs. This step may take a while to finish.

To get help documentation about this module, run:

 perl lncRNApipe --h cpc

=item 6. RNAfold (lncRNApipe Option --rna or --rnafold)

Here, lncRNApipe pipeline invokes RNAfold program to calculate minimum free energy structure of the predicted non-coding RNA. This step may take a while to finish.
                                                                                                                                                                
To get help documentation about this module, run:
    
 perl lncRNApipe --h rna

=item 7. Infernal (lncRNApipe Option --inf or --infernal)

In the final step, lncRNApipe pipeline invokes cmscan from Infernal package to search Rfam databases for sturcture and sequence similarities to annotate the putative lncRNAs .
                                                                                                                                                                
To get help documentation about this module, run:
    
 perl lncRNApipe --h inf

=back

=head1 CAVEATS

The pipeline script uses a lot of inherent GNU core utils and has been only tested in BASH shell. Please use absolute full PATH names.

=over 5

=item * Instead of using:

  lncRNApipe --run ./lncRNApipe_output ...

 Use: 

  lncRNApipe --run /data/lncRNApipe_output ...

=back

=head1 OUTPUT

The final output is a putative list of annotated lncRNAs in GTF format that can be directly uploaded as custom tracks to Genome Browsers, such as UCSC etc... and predicted secondary structures by RNAfold. 

=over 5

=item * Output from "cuffcompare" is stored in the directory called "cuffcompare". lncRNApipe uses "lncRNApipe_cuffcmp.tracking" file from this folder for the next module.

=item * Output from "categorize_ncRNAs.pl" is stored in the directory called "categorize_ncRNAs". lncRNApipe uses files ending with suffix ".putative.class.lncRNAs.gtf" for the next module.

=item * Output from "get_unique_features.pl" output is stored in the directory called "get_unique_features". lncRNApipe uses files ending with suffix ".putative.class.lncRNAs.unique.gtf" for the next module

=item * Output from "fetch_seq_from_ucsc.pl" is stored in the directory called "fetch_seq_from_ucsc". lncRNApipe uses these FASTA seqeunce files to run CPC, RNAfold and Infernal modules. The file ending in suffix ".noncoding.FASTA" is used for the next modules.

=item * CPC predictions are stored in a directory called "CPC".

=item * The final result files are stored in a directory called "lncRNApipe.final" and have a suffix ".final.gtf".

=item * The RNAfold plots are stored in directories inside "lncRNApipe.final" directory ending in suffix ".final.RNAfold".

=item * The result files from "cmscan" (Infernal) are stored in directories inside "lncRNApipe.final" directory ending with suffix ".infernal". The file ending in suffix ".allInfernalHits.gtf" contains all RNA hmm and cm model matches from "cmscan" (Infernal) and the best hit is applied to file ending in suffix ".final.gtf" based on the sequence coverage ("--coverage-infernal") filter requested.

=back

The pipeline has the ability to run each of these individual modules if any of them fail during the initial run.

=head1 EXAMPLE

=head2 Installation:

=over 5

=item * This script will try its best to setup the pipeline on the suitable architecture.

 perl lncRNApipe --setup

=back

=head2 Start the Pipeline:

=over 5

=item * Run the complete pipeline:

 perl lncRNApipe --run /data/lncRNApipe/ --cuff '-r /data/projects/mm10/macrophages/rna-seq/cuffcompare.m0.m1.m2_vs_ensembl_ref_known/refSeq_UCSCKnownGenes_Ensemble.gtf -s /data/ref_fasta/mm10_UCSC/whole_genome.fa /data/projects/mm10/macrophages/rna-seq/m0/cufflinks/transcripts.gtf /data/projects/mm10/macrophages/rna-seq/m1/cufflinks/transcripts.gtf /data/projects/mm10/macrophages/rna-seq/m2/cufflinks/transcripts.gtf' -cat '-clean -overlap 80 -len 200 -max-len 10000 -min-exons 1 -fpkm 2 -antisense -sample-names "M0,M1,M2"' --get '-sf /data/projects/mm10/macrophages/ncrna/known_ref_lncRNAs/mm10_ncRNA.bed' --fetch --cpc --inf

=item * Run categorize_ncRNAs module of the pipeline:
    
 perl lncRNApipe --run /data/lncRNApipe/ --cat '-clean -overlap 80 -len 200 -max-len 10000 -min-exons 1 -fpkm 2 -annotation /data/projects/mm10/macrophages/ncrna/2lncRNA/allGenes.txt -antisense -sample-names "M0,M1,M2" /data/projects/mm10/macrophages/rna-seq/m0/cufflinks/transcripts.gtf /data/projects/mm10/macrophages/rna-seq/m1/cufflinks/transcripts.gtf /data/projects/mm10/macrophages/rna-seq/m2/cufflinks/transcripts.gtf'

=item * Run get_unique_features module of the pipeline to extract novel lncRNAs:
    
 perl lncRNApipe -run /data/lncRNApipe --get '-sf /data/projects/mm10/macrophages/ncrna/known_ref_lncRNAs/mm10_ncRNA.bed'

=item * Run get_unique_features module of the pipeline to extract known lncRNAs:
    
 perl lncRNApipe --run /data/lncRNApipe --get '--known -sf /data/projects/mm10/macrophages/ncrna/known_ref_lncRNAs/mm10_ncRNA.bed'

=item * Run fetch_unique_features module of the pipeline:

 perl lncRNApipe --run /data/lncRNApipe --fetch '-db mm10'

=item * Run CPC module of the pipeline:

 perl lncRNApipe --run /data/lncRNApipe --cpc

=item * Run RNAfold module of the pipeline:

 perl lncRNApipe --run /data/lncRNApipe --rna

=item * Run Infernal module of the pipeline:

 perl lncRNApipe --run /data/lncRNApipe --inf

=back

=head1 OPTIONS

lncRNApipe takes the following arguments:

=over 4

=item --v or --version (Optional)
    
Displays version information and quits.

=item --d or --no-update-chk (Optional)
    
Disable software update check.

=item --h or --help (Optional)
    
Displays this helpful message.

=item --setup or --setup
    
Try to setup the pipeline with all its dependencies.

=item --run or --run

Run the lncRNApipe pipeline with output directory supplied as option.
Ex: perl lncRNApipe -run /home/konganti/lncRNApipe_output ...

=item --cuff or --cuffcompare

Run the Cuffcompare module with supplied options.

=item --cat or --cat-ncRNAs

Run the categorize_ncRNAs.pl module with supplied options.

=item --get or --get-uq-feat

Run the get_unique_features.pl module with supplied options.

=item --fetch or --fetch-seq

Run the fetch_seq_from_ucsc.pl module with supplied options.
Using --fetch without any options will try and use local reference FASTA
mentioned with --cuffcompare option.

=item --skip-get (--fetch-seq option required)

Run the fetch_seq_from_ucsc.pl module without running get_unique_features.pl module.
In some cases, you may not have a file of known non-coding RNAs. In such cases,
you can skip get_unique_features.pl module, and directly fetch FASTA sequences for
features from --categorize_ncRNAs.pl module.

=item --cpc or --cpc

Run the CPC module. No options are required.

=item --skip-cpc

Skip runnning CPC altogether. Use this option, if you think CPC is flagging a lot of
your transcripts as "coding" and instead rely only on Infernal search results.

=item --skip-blastall-core

Skip runnning core CPC process once you know you have output from CPC.
This option can be used when lncRNApipe fails for some reason after blastall
within CPC has finished running but unable to continue forward.

=item --rna or --rnafold

Run RNAfold module of the pipeline. No options required

=item --skip-rnafold-core

Skip runnning core RNAfold process once you know you have output from RNAfold.
This option can be used when you want lncRNApipe to generate plots based on
output from RNAfold which should have completed successfully.

=item --inf or --infernal

Run Infernal module of the pipeline. No options required

=item --skip-cmscan-core

Skip runnning cmscan process once you know you have output from a successful 
cmscan run. This option can be used when you want lncRNApipe to extract 
annotation from infernal and attempt to annotate putative novel lncRNAs.

=item --cov-inf or --coverage-infernal

Only annotate final putative lncRNAs which have cmscan match over this much percentage of sequence.
Default: 10

=item --cpu

lncRNApipe will attempt to run multiple processes on this may CPU cores.
Mentioning number of CPUs equal to or more than number of assembled 
transcript files is strongly encouraged.

=item --send-run-report

lncRNApipe will send your log file. For this to work, you need to save the output of the pipeline to
a log file. Ex:

     perl lncRNApipe --run /data/lncRNApipe/ \
                     --cuffcmp '-r annotation.gtf -s genome.fa transcripts1.gtf transcripts2.gtf' \
                     --cat-ncRNAs '-sample-names "M1,M2" -ov 80 -fpkm 2 -len 200 -max-len 10000 -min-exons 1 -antisense" \
                     --fetch-seq '-db mm10' \
                     --cpc \
                     --rna \
                     --inf  &> lncRNAPipe.run.log 

and then, if the run fails or if you wish to discuss other issues. Use:
 
     perl lncRNApipe --send-run-report '-email your_email@gmail.com -log lncRNApipe.run.log -m "I cannot run lncRNApipe because lorem ipsum ...."'

=back

=head1 AUTHOR

Kranti Konganti, E<lt>konganti@tamu.eduE<gt>.

=head1 COPYRIGHT

This program is distributed under the Artistic License.

=head1 DATE

May-08-2015

=cut
