#!/bin/bash

# This script will automaically execute erveryday @3:59 AM
# edit file (sudo) /etc/crontab to change the time
# 59 3 * * * root bash /home/wangming/bin/update_FTP_files.sh >> /data/yulab_ftp/.log/update.stats

# change the permission and owner:group of files in /data/yulab_ftp/Share
#
# Files:
# permission: rw-r--r--
# ownership:  wangming:yulab
#
# Folders:
# permission: rwxr-xr-x
# ownership:  wangming:yulab
# 2016-10-26

function update_file {
    FT=(`stat -c "%A %U %G" "$1"`) # drwxr-xr-x wangming yulab
    FRIGHT=${FT[0]}
    FOWNER=${FT[1]}
    FGROUP=${FT[2]}
    ## Target
    TRIGHT="rwxr-xr-x"
    T1="pass" ; T2="pass" ; T3="pass"
    if [[ $FRIGHT == d* ]] ; then 
    # Directory
        [[ $FRIGHT != drwxr-x--- ]] && chmod 750 $1 && T1="${FRIGHT}"
        [[ $FOWNER != yuftp ]] && chown yuftp $1 && T2="${FOWNER}"
        [[ $FGROUP != yulab ]] && chgrp yulab $1 && T3="${FGROUP}"
        echo -e "D: $T1\t$T2\t$T3\t$1"
    elif [[ $FRIGHT == -* ]] ; then
    # Regular file
        [[ $FRIGHT != -rw-r----- ]] && chmod 640 $1 && T1="${FRIGHT}"
        [[ $FOWNER != yuftp ]] && chown yuftp $1 && T2="${FOWNER}"
        [[ $FGROUP != yulab ]] && chgrp yulab $1 && T3="${FGROUP}"
        echo -e "F: $T1\t$T2\t$T3\t$1"
    elif [[ $FRIGHT == l* ]] ; then
    # Symbolic link
        chown 777 $1
        echo -e "L: $T1\t$T2\t$T3\t$1"
    else
        echo -e "X: Not recognized\t$1"
    fi
}

############
## config ##
############
BASE_DIR="/data/yulab_ftp"
#BASE_DIR="/home/wangming/work/test"
SHARE_DIR=${BASE_DIR}/Share
#SHARE_DIR=${BASE_DIR}/02.Group_meetings

##
DT_Y=`date "+%Y"`
DT_M=`date "+%m"`
DT_D=`date "+%d"`
DT_T=`date "+%H%M%S"`
LOG_DIR="${BASE_DIR}/.log/${DT_Y}/${DT_M}"
LOG_FILE="${LOG_DIR}/${DT_Y}${DT_M}${DT_D}_${DT_T}.log"
LOG_TEMP="${LOG_DIR}/${DT_Y}${DT_M}${DT_D}_${DT_T}.temp"
LOG_ERR="${LOG_DIR}/${DT_Y}${DT_M}${DT_D}_${DT_T}.err"
[[ ! -d $LOG_DIR ]] && mkdir -p $LOG_DIR
## search - 1st
find $SHARE_DIR > $LOG_TEMP 2> $LOG_ERR
## update files
while read LINE
do
    # consider blanks in fliename
    echo "$LINE" | rename -e 's/ |\?/./g' 
    update_file "$LINE" >> $LOG_FILE
done < $LOG_TEMP

## update folders that failed to access
for i in $(seq 1 4) # go further 4-level in that directory
do
    DT_F=`date "+%Y-%m-%d %H:%M:%S"`
    ERR_LINE=`wc -l $LOG_ERR | cut -d" " -f 1`
    [ $ERR_LINE == 0 ] && echo "## Finish @ $DT_F" && break 
    DIR_N=()
    while read DL
    do
        DN=($DL)
        DT=${DN%\'}
        DT=${DT#\`}
        DIR_N+=($DT)
    done < $LOG_ERR
    EXT_DD=$(IFS=" "; echo "${DIR_N[*]}")
    find $EXT_DO > ${LOG_TEMP}.$i 2> $LOG_ERR
    while read LINE
    do
        update_file $LINE
    done < ${LOG_TEMP}.$i
    rm ${LOG_TEMP}.$i $LOG_ERR
done
rm $LOG_TEMP

## EOF
