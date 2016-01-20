#!/bin/bash

#################
#    Compare    #
#################

### input a name of path
perl bin/anno_files_sRNAs.v3.pl -i db/Kroger_2013_STyph_280sRNA_list.txt data1/01_EEP

### input a list of paths 
perl bin/anno_files_sRNAs.v3.pl -i db/Kroger_2013_STyph_280sRNA_list.txt tag.list1

### input a list of paths, (alternative directory structure)
perl bin/anno_files_sRNAs.v3.pl -i db/Kroger_2013_STyph_280sRNA_list.txt tag.list2

### support different type of files, [sRNA|UTR|IM]
perl bin/anno_files_sRNAs.v3.pl -i db/Kroger_2013_STyph_280sRNA_list.txt -t UTR tag.list2

### change the chr_name of both db and query to the same
perl bin/anno_files_sRNAs.v3.pl -i db/Kroger_2013_STyph_280sRNA_list.txt -t UTR -n chr tag.list2

#################
# Visualization #
#################




