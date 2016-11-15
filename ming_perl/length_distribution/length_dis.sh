#!/bin/bash

# calculate the length distribution of reads in FASTQ format

function usage () {
cat << EOF
Usage: length_dis.sh <fastQ>

Calculate the length distribution of reads in FASTQ format.
EOF
exit 0
}

if [ $# -lt 1 ] ; then 
    usage
fi

INFILE=${1}

if [ -f ${INFILE} ] ; then
    cat ${INFILE} | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths){print l, lengths[l]}}'
else
    echo "${INFILE} file not exists"
    exit 1
fi



