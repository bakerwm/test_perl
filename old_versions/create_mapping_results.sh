From Clean fasta to soap & match_genome

1. Soap >> (Bi server, 1 cpu for each job)
	(1) 45SE & 81SE: -M 0  -r 2  -v 0   -p 1
example: soap -a clean.fa  -D NC_000962.fna.index -o Rv_45SE.soap  -u Rv_45SE.unmap  -M 0  -r 2  -v 0   -p 1
	
	(2) 90PE:                -m2   -r 2  -v 2   -p 1   -x 1000  -s 40  -l 32  -o out.soapPE  -2 out.soapSingle
example: soap -a Rv_1.fq  -b Rv_2.fq  -D NC_000962.fna.index -o Rv_200PE.soap  -2 Rv_200PE.soapSingle  -m2   -r 2  -v 2   -p 1   -x 1000  -s 40  -l 32 

2. Trans Soap:
	(1) match_genome.txt [7-line tab-delimited file]: <tag_id> <genome_id> <begin> <end> <strand> <sequence> <exp>
a. SE soap:
"soap2txt_v1.pl"
usage: soap2txt_v1.pl  clean.fa  *.soap
	
b. PE soap:
"soapPE2SE.pl"
usage: 

merge.pl Rv_L1.PESoap 1>Rv_L1_genome.txt 2>length_vs_number_Rv_L1_genome.txt

	(2) coverage [3-line tab-delimited file]: <genome_id> <position> <coverage>
"soap.coverage"
usage: soap.coverage  -cvg  -i  *.soap  -refsingle  H37Rv.fna  -o  result.txt  -depthsingle  H37Rv_45SE.cov  
"CovTrans.pl"
Usage: perl  CovTrans.pl  in.cov   result.out  Genome_fna  out.coverage
Note:
1. Infile is the result from soap.coverage -depthsingle
2. Length is the result from soap.coverage -o
3. The genome fasta file used for soap
4. The result file
	
3. Get seed sequence:
	(1) get seed sequence from coverage file.
"get_lncRNA_v2.pl"
usage: perl  get_lncRNA_v2.pl  -i coverage  -s strand  -c cutoff  -o outfile
Options:
    -i <file>       : coverage file, 3-line tab-delimited file in the following format:
                      <Genome_ID> <Position> <Coverage>
                      Position should be continuous from 1 to end.
    -s [+/-]        : The strand.
    -c [0-1]        : The cutoff [0-1] to determin the edge. Default 0.5 .
    -o <file>       : The result file.

	(2) Format transformation.\
"getPosition.pl"
Usage: perl /home/wangming/work/bin/get_sRNA/getPosition.pl  infile gff  strain_id  > outfile
Note:
    -1      : The input file in the following format:
              <ID> <Strand> <Begin:Cov> <Max:Cov> <End:Cov>
    -2      : The GFF file. GFF3 recommend.
    -3      : The Strain name, help to recognize the GFF file.
Example: perl  getPosition.pl  H37Rv_45SE.temp  NC_000962.gff H37Rv  >H37Rv_45SE_lncRNA.txt

	(3) Calculate rpkm in each library.
"?"
...


	(4) Merge sRNA files from different libraries.
"merge4all.pl"
Usage: perl merge4all.pl -n 4 -i in_1,in_2,in_3,in_4 -o outfile
Options:
    -n <Integer>    :the number of input files
    -i <file>       :input files
    -o <file>       :output file
    -h              :show this help

Note: For multiple input files,
1. Filename like this, "45SE_H37Rv_*". eg: 45SE_H37Rv.txt
2. Multiple input filenames join by comma ",". eg: "infile_1,infile_2"

Description:
1. Two sequences have an overlap longer than 1/2 of one of the 
   sequence or longer than 20 nt will merge into a new one.
   (match_length > 1/2 of input or 20 nt)
2. The expression level between sequences should not more than 
   100-fold.
   (Max_exp/exp <= 100)
   
	(5) RNA classification.
"sort2candi_v1.pl"
Usage: perl sort2candi_v1.pl  infile

Note: the rules to select sRNA
1. at IGR or AS.
2. >=100 bp from 5' ORF and >=60 bp from 3' ORF.


==========================================================================================
End of file
==========================================================================================

>> Used command lines.

#command lines.
From  rpkm.txt to merge lncRNA.txt.
	get_work_shell.pl
	
	[getPosition.pl, get_sRNA_sh.pl  sort2candi_v1.pl  sort2temp.pl]
	
From FastA files to rpkm.txt.


From match_genome to mapping.
(1)
	Cal_sRNA_end.pl
	getRscript.pl
	Rscript  *.R
	
(2)
	Cal_sRNA_end.pl
	match_v2.pl
	getRscript_v2.pl
	Rscript *.R

perl  ~/work/bin/merge.pl  200PE_rpkm.txt  140PE_rpkm.txt  81SE_rpkm.txt  45SE_rpkm.txt   >merge.txt
perl ~/work/bin/pick_merge_sRNA.pl   merge.txt  >merge.temp

perl  ~/work/bin/getPosition.pl   NC_000962.gff   merge.temp  >merge_lncRNA.txt

perl ~/work/bin/sort2candi_v3.pl     merge_lncRNA.txt 

# check reads mapping coverage on each mRNA.

Declaration: 
1. Prepare samples.
2. Change file name.
3. Confirm parameters.

Step 1. Clean fasta file. 
1. SOAP or BWA or Bowtie
2. match_genome.txt
3. Poscoverage.txt
4. Match genome/mRNA statment.

Step 2. Get sRNAs.
1. Statistic poscoverage, 5', 3' and read.
2. Merge different libraries. Twice.
3. Get position.
4. Find sRNAs.

Process:
1. Begin with:
	a. match_genome.txt files,
	b. Genome annotation files: GFF, fna,
2. make dir for match_genome and annotation files.
	$ mkdir  H37Rv_match_genome/	 database/		# change file name to "match_genome.45SE.txt";
	$ perl run_statPoscov.pl     H37Rv_match_genome/
	"45SE"_3end.mapping.txt

3. Merge sRNAs from different libraries.
	$ perl  merge4all  -n 4	-i  files,files  -o outfile
	$ perl  sort2temp.pl  outfile 
	$ perl  getPosition		GFF  outfile.temp
	$ perl  sort2candi_v3.pl  outfile.txt

Step 3. Statistic sRNAs mapping tRNA/rRNA/pub_RNA...
	a. draw input file coverage map.
		$ mkdir  mRNA_maps
		$ ln -s  genome_mapping.txt ./
		$ perl  stat3end.pl
		$ perl  stat5end.pl
		$ perl  statPoscov.pl
		$ perl  creat_title_dir.pl
		$ perl  getRscript.pl
		$ Rscript *R
		
	b. draw sRNAs mapping input files.
		$ perl  mkdir mRNA_match_sRNA
		$ perl  match.pl  sRNAs  inputfile
		$ perl  stat3end.pl
		$ perl  stat5end.pl
		$ perl  statPoscvo.pl
		$ perl  creat_title_dir.pl
		$ perl  getRscript_v2.pl  infile   total_line  total_seed
		$ Rscript   input2pdf.R
rm di
Step 4. Summary
1. Match genome/mRNA/rRNA/tRNA/IGR/AS/ statment.
2. sRNA output.
3. Draw figures: match mRNA/rRNA/tRNA...
