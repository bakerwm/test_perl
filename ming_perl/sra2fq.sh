
# convert SRA to fastq files

if [ $# != 1 ]; then
    echo 'Usage: bash sra2fq.sh <sra.dir>'
    echo 'Save fastq files to the same directory as *.sra'
    exit
fi

SRA_DIR=`readlink -f $1`

for i in ${SRA_DIR}/*.sra
do
    if [ -f $i ]; then
        echo $i
        fastq-dump $i
    fi
done

