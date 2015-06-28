set -ue

## The following procedure described in:
# http://khmer-protocols.readthedocs.org/en/v0.8.4/mrnaseq/index.html

## The overall process
## 1. Trimmomatic - cut adapters
## 2. khmer - digitally normalize reads
## 3. Trinity - Assemble reads
## 4. eel-pond - annotate assembly
## 5. RSEM + EBSeq - Differential expression

## required softwares
# 1. screed
# 2. khmer
# 3. Trimmomatic
# 4. libgtextutils and fastx
# 5. Trinity
# 6. bowtie
# 7. samtools
# 8. BLAST
# 9. rsem

## data & tools
data=/home/wangming/Identify_sRNA/khmer-protocols/data
EcoliFaa=/home/wangming/work/database/E.coli/NC_000913.faa
trimDir=/home/wangming/Documents/Trimmomatic-0.32
khmerDir=/usr/bin
epDir=/home/wangming/Identify_sRNA/khmer-protocols/software/eel-pond
rsemDir=/home/wangming/Identify_sRNA/khmer-protocols/software/RSEM-1.2.20
trinityDir=/home/wangming/Identify_sRNA/khmer-protocols/software/trinityrnaseq-2.0.6
#
trimJar=${trimDir}/trimmomatic-0.32.jar
trimAD=${trimDir}/adapters/TruSeq2-PE.fa
interLvRd=${khmerDir}/interleave-reads.py
splitPE=${khmerDir}/split-paired-reads.py
doPart=${khmerDir}/do-partition.py
normByMedian=${khmerDir}/normalize-by-median.py
filterAbund=${khmerDir}/filter-abund.py
extractPE=${khmerDir}/extract-paired-reads.py
Trinity=${trinityDir}/Trinity
renaPart=${epDir}/rename-with-partitions.py
makeUnBH=${epDir}/make-uni-best-hits.py
makeReBH=${epDir}/make-reciprocal-best-hits.py
makeNmdb=${epDir}/make-namedb.py
AnnoSeq=${epDir}/annotate-seqs.py
makeTrGeMpFl=${epDir}/make-transcript-to-gene-map-file.py
extAnCh=${epDir}/extract-and-annotate-changed.py
plotEp=${epDir}/plot-expression.py
rsemPrRf=${rsemDir}/rsem-prepare-reference
rsemCaEp=${rsemDir}/rsem-calculate-expression
rsemGeDaMx=${rsemDir}/rsem-generate-data-matrix
rsemRnEb=${rsemDir}/rsem-run-ebseq

##################
## 1/8. trim data
##################

## make wkdir
mkdir -p 1.trim
cd 1.trim

## data
for i in $data/*.fastq
do
    echo "$(basename $i)"
    head -40000 $i | gzip > $(basename $i)\.gz
#    gzip -c $i > $(basename $i).gz
done

### trim adaptors
for i in *_1.fastq.gz
do
    R1=../$i
    R2=${R1//_1/_2}
    NAME="$(basename $i _1.fastq.gz)"    
    # make a temp directory
    mkdir -p trim
    cd trim    
    # run trimmomatic
    java -jar $trimJar PE $R1 $R2 s1_pe s1_se s2_pe s2_se ILLUMINACLIP:$trimAD\:2:30:10 
    # change the id of fastq
    awk '{if(NR%4==1) $1=$1"/1"; print;}' s1_pe > s1_pe.id
    awk '{if(NR%4==1) $1=$1"/2"; print;}' s2_pe > s2_pe.id
    # interleave the remaining paired-end files
    $interLvRd s1_pe.id s2_pe.id > ../$NAME.pe.fq
    gzip -f ../$NAME.pe.fq
    # combine the single-ended files
    cat s1_se s2_se | gzip -fc > ../$NAME.se.fq.gz
    # go back up to the working directory and remove the temp directory
    cd ..
    rm -r trim
    # make it hard to delete the files you just created
    chmod u-w *.pe.fq.gz *.se.fq.gz
    # 2. quality trim (skipped)
    cp $NAME.pe.fq.gz $NAME.pe.qc.fq.gz
    cp $NAME.se.fq.gz $NAME.se.qc.fq.gz
    # 3. Extracting paired ends from files
done

#chmod u-w *.qc.fq.gz

# 4. remove other files
rm -fr *.fastq.gz *.pe.fq.gz *.se.fq.gz

cd ..

############################
# 2/8. normalizing data
############################

## make wkdir
mkdir -p 2.norm
cd 2.norm

## data
ln -fs ../1.trim/*.qc.fq.gz .

## run digital normalization 
$normByMedian -q -p -k 20 -C 20 -N 4 -x 2e9 -s table.ct *.pe.qc.fq.gz
$normByMedian -q -C 20 -l table.ct -s table.ct *.se.qc.fq.gz

## trim off likely erroneous k-mers
$filterAbund -C 2 -V table.ct *.keep

## rename files
# PE, break out the orphaned and still-paired reads:
for i in *.pe.*.abundfilt
do
    $extractPE -f $i
done

# SE, combine the orphaned reads into a single file:
for i in *.se.*.abundfilt
do
    pe_orphans=$(basename $i .se.qc.fq.gz.keep.abundfilt).pe.qc.fq.gz.keep.abundfilt.se
    newfile=$(basename $i .se.qc.fq.gz.keep.abundfilt).se.qc.keep.abundfilt.fq.gz
    cat $i $pe_orphans | gzip -c > $newfile
done

# rename the remaining PE reads & compress those files
for i in *.abundfilt.pe
do 
    newfile=$(basename $i .fq.gz.keep.abundfilt.pe).keep.abundfilt.fq
    mv $i $newfile
    gzip $newfile
done

rm -fr *.se.*.abundfilt *.pe.*.abundfilt *.pe.*.abundfilt.se *.keep *.abundfilt *.qc.fq.gz *.ct

cd ..

########################
# 3/8. assembly
########################

## make wkdir
mkdir 3.assembly
cd 3.assembly

## data
ln -fs ../2.norm/*.fq.gz .

# build the files to assemble
for i in *.pe.qc.keep.abundfilt.fq.gz
do
    $splitPE $i
done

cat *.1 > left.fq
cat *.2 > right.fq

gunzip -c *.se.qc.keep.abundfilt.fq.gz >> left.fq

## Run trinity
$Trinity --seqType fq --max_memory 10G --left left.fq --right right.fq --CPU 10

cd ..

#########################
# 5/8. build transcripts
#########################

## make wkdir
mkdir 5.build_trans
cd 5.build_trans

## data
gzip -c ../3.assembly/trinity_out_dir/Trinity.fasta > trinity-nematostella-raw.fa.gz

## Run khmer partitioning
$doPart -x 2e9 -N 4 --threads 4 nema trinity-nematostella-raw.fa.gz

python $renaPart nema trinity-nematostella-raw.fa.gz.part
mv trinity-nematostella-raw.fa.gz.part.renamed.fasta.gz trinity-nematostella.renamed.fa.gz

cd ..

###############################
# 6/8. annotating transcripts
###############################

## make wkdir
mkdir -p 6.annotating_transcripts
cd 6.annotating_transcripts

## data
cp -r ../5.build_trans/trinity-nematostella.renamed.fa.gz .

## unzip transcript seq
gunzip trinity-nematostella.renamed.fa.gz

## Download reference RefSeq
ln -fs $EcoliFaa mouse.protein.faa

## Build blast db
formatdb -i mouse.protein.faa -o T -p T
formatdb -i trinity-nematostella.renamed.fa -o T -p F

## run blast
blastall -i trinity-nematostella.renamed.fa -d mouse.protein.faa -e 1e-3 -p blastx -o nema.x.mouse -a 8 -v 4 -b 4
blastall -i mouse.protein.faa -d trinity-nematostella.renamed.fa -e 1e-3 -p tblastn -o mouse.x.nema -a 8 -v 4 -b 4

## reciprocal best hits
python $makeUnBH nema.x.mouse nema.x.mouse.homol
python $makeReBH nema.x.mouse mouse.x.nema nema.x.mouse.ortho

## prepare some mouse info
python $makeNmdb mouse.protein.faa mouse.namedb 
python -m screed.fadbm mouse.protein.faa

## annotate sequences
python $AnnoSeq trinity-nematostella.renamed.fa nema.x.mouse.ortho nema.x.mouse.homol

## rename output file
cp trinity-nematostella.renamed.fa.annot nematostella.fa

cd ..

################################
# 7/8. expression analysis RSEM
################################

## make wkdir
mkdir -p 7.exp_rsem
cd 7.exp_rsem

## data
ln -fs ../6.annotating_transcripts/nematostella.fa .

## make a tran_to_gene_map file
python $makeTrGeMpFl nematostella.fa nematostella.fa.tr_to_genes

## ask RSEM for reference
$rsemPrRf --bowtie -q --transcript-to-gene-map nematostella.fa.tr_to_genes nematostella.fa nema

## find and list reads
ln -fs ../1.trim/*.pe.qc.fq.gz .

## make a list of data files
ls *.pe.qc.fq.gz > list.txt

## Run RSEM    
n=1
for filename in $(cat list.txt)
do
    echo mapping $filename
    gunzip -c $filename > ${n}.fq
    split-paired-reads.py ${n}.fq
    $rsemCaEp --paired-end ${n}.fq.1 ${n}.fq.2 nema -p 4 ${n}.fq
    rm ${n}.fq ${n}.fq.[12] ${n}.fq.transcript.bam ${n}.fq.transcript.sorted.bam
    n=$(($n + 1))
done

## gather results
$rsemGeDaMx [1-6].fq.genes.results > a-vs-b.matrix

## other results
$rsemGeDaMx 1.fq.genes.results 3.fq.genes.results > results.matrix

cd ..

##########################
# 8/8. differential EBSeq
##########################

## make wkdir
mkdir -p 8.ebseq
cd 8.ebseq

## data
cp ../7.exp_rsem/a-vs-b.matrix .
ln -fs ../7.exp_rsem/nematostella.fa .

## differential expression with EBSeq

## input data
# nematostella.fa, [0-9].fq.genes.results (a-vs-b.matrix)

## run EBSeq
$rsemRnEb a-vs-b.matrix 3,3 a-vs-b.changed

## extract diff_exp_genes
python $extAnCh a-vs-b.changed nematostella.fa a-vs-b.changed.csv

## visualize diff_exp_genes
python $plotEp a-vs-b.matrix 3,3 a-vs-b.changed.csv

## output a *.png file.
cd ..

## __END of Script__ ##
