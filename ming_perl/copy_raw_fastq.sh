#!/bin/bash 

## copy and rename the raw fastq files from "seq_data" 
## to your working dir
##
## rename the file according to the file: readme.txt in the same dir
##
## <SRR000001> <name>
##
## rename SRR000001.fastq to name.fastq in you working dir

if [ $# != 1 ]; then
    echo 'Usage: copy_raw_fastq.sh <fq.dir>'
    echo 'copy and rename the fastq file according to the readme.txt file'
    exit
fi

## copy and rename the data to this dir
DT_DIR=$1
#DT_DIR="$HOME/seq_data/CLIP/iCLIP/hnRNP_C"
README="${DT_DIR}/readme.txt"

if [ -f $README ]; then
    ln -fs ${README} .
else 
    echo 'Cannot find readme.txt in: $DT_DIR'
fi

##
echo "[Copy and rename FASTQ files]"

while read LINE
do
    a=( $LINE )
    OLD=${DT_DIR}/${a[0]}.fastq
    NEW=${a[1]}.fastq
    if [ -f ${OLD} ]; then
        echo -e "create ${NEW} from ${OLD}"
        ln -fs ${OLD} ${NEW}
    fi
done < ${README}

## EOF
