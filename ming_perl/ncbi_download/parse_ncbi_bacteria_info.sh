#!/bin/bash
#
# Download bacteria fna files
#
#
. ${HOME}/software/src/bash/utilities.sh
baseurl="ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/"

function list_ftp {
    url=$1
    curl --silent -l ${url}/
}

function parse_fna {
    url=$1
    ids=($(list_ftp ${url} | xargs))
    for id in ${ids[@]} ; do 
        [[ ${id} == *.fna ]] && echo ${id%.fna}
    done
}

function parse_bacid {
    bacid=$1
    bacurl="${baseurl}/Bacteria/${bacid}/"
    fnas=($(parse_fna ${bacurl} | sort | xargs))
    echo "${bacid} ${fnas[@]}"
}

function parse_bac {
    url="${baseurl}/Bacteria/"
    curl --silent -l ${url} | grep "_uid"
}

function download_file {
    url=$1
    outdir=$2
    prefix=$(basename ${url})
    [[ ! -d ${outdir} ]] && mkdir -p ${outdir}
    curl --silent -o ${outdir}/${prefix} ${url}
}

## main
[[ $# -lt 2 ]] && echo "Usage: $0 <out.dir> <pattern> <download_files:0|1>" && exit 0
outdir=$1
patt=$2
download=$3
outinfo="${outdir}/bacteria_fna.list"
[[ ! -d ${outdir} ]] && mkdir -p ${outdir}
[[ -s ${outinfo} ]] && echo "${outinfo} - file is not empty, exiting ..." && exit 1 
ids=($(parse_bac | sort | grep "${patt}" | xargs))
[[ ${#ids[@]} == 0 ]] && echo "[${patt}] - not hit bacteria genomes, see: ${baseurl}/Bacteria/" && exit 1
[[ ${#ids[@]} -gt 1000 ]] && echo "[${patt}] - too much hits" && exit 1
p=10
n=0
pids=()
#echo ${ids[@]}
for id in ${ids[@]} ; do 
    parse_bacid ${id} >> ${outinfo}&
    pids+=("$!")
    n=$((n + 1))
    echo2 "${id} [${n}/${#ids[@]}]"
    if [[ $[${n} % ${p}] == 0 ]] ; then
        wait ${pids[@]}
        pids=()
    fi
done
wait ${pids[@]} # wait more proc

## check, whether to download files
[[ -z ${download} ]] && download=0
if [[ ${download} == 1 ]] ; then
    while read line ; do 
        fs=(${line})
        name=${fs[0]}
        unset fs[0]
        for f in ${fs[@]} ; do 
            echo2 "${f} downloading ..."
            furl="${baseurl}/Bacteria/${name}/${f}.fna"
            wget -q -P ${outdir} ${furl}
        done
    done < ${outinfo}
fi


## EOF
