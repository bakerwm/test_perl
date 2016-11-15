#!/bin/bash

#####################################################################
# This script is designed to calculate the kmer frequency for a     #
# input fasta file                                                  #
# output with the following format                                  #
# <mer> <frequency>                                                 #
#                                                                   #
# Wang Ming 2016-11-15                                              #
#####################################################################

## parameters
## <mer_length> <file>
if [ $# != 2 ] ; then
    echo "Usage: $0 <mer_length> <file>"
    exit
fi
MER=$1
FILE=$2
[[ ! -f $FILE ]] && (echo "File not exists - $FILE" ; exit)
if [[ $FILE != *.f[aq]* ]] ; then echo -e "$FILE - File not supported; \nOnly accept: fa|fq|fasta|fastq, Exiting ..."; exit ; fi

## required tools
## faSize, jellyfish
function checkBin {
    which "$1" &> /dev/null
    if [[ $? != 0 ]] ; then echo "Cannot find software ${1}, Exiting..." ; exit ; fi
}
checkBin "md5sum"
checkBin "faSize"
checkBin "jellyfish"

## estimate the size of fa
FA_SIZE=`faSize $FILE | head -n 1 | awk '{print $1}'`
RND_NUM=`md5sum $FILE | cut -d" " -f 1`
TEMP_JF="temp.${RND_NUM}.jf"
jellyfish count -m $MER -s $FA_SIZE -t 8 -o $TEMP_JF $FILE
jellyfish dump -c -t -L 0 $TEMP_JF | awk -v OFS="\t" '{a[$1]=$2; sum+=$2}END{for (i in a) {print i,a[i],(a[i]/sum)}}'
rm $TEMP_JF

## EOF
