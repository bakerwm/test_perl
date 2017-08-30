#!/bin/bash

# This script will automaically execute erveryday @3:59 AM
# edit file (sudo) /etc/crontab to change the time
# 59 3 * * * root bash /home/wangming/bin/update_FTP_files.sh >> /data/yulab_ftp/.log/update.stats

# change the permission /data/yulab_ftp/Share (rwxrwx---)
#
# Files:
# permission: rw-r-----
# ownership:  wangming:yulab
#
# Folders:
# permission: rwxrwx---
# ownership:  wangming:yulab
#
# Backup (hourly)
# 1. create a backup archive
# 2. update only new files, and create a backup
#
# 2016-12-02 version 0.3
#   - add backup function
#   - hourly update (*:59)
# 2016-10-26 version 0.2
#   - rewrite update function
#   - daily update (4:00 am)

############
## config ##
############
BASE_DIR="/data/yulab_ftp/"
#SHARE_DIR=${BASE_DIR}test
SHARE_DIR=${BASE_DIR}Share
BACKUP_DIR="${BASE_DIR}Archive/Share_backup"
#SHARE_DIR=${BASE_DIR}02.Group_meetings
##
[ ! -d $BACKUP_DIR ] && mkdir -p $BACKUP_DIR
## update $BACKUP_DIR permission
chmod 750 $BACKUP_DIR
chown yuftp:yulab $BACKUP_DIR

## functions
## update files
function update { 
    ## input: <infile> <0|1>
    ## update owner or not
    infile=$1
    ## permission of folder 0=770, 1=750
    if [ x$2 != x ] ; then
	dir_mode=$2
    else
	dir_mode=0
    fi
    ## owner of file/dir
    if [ x$3 != x ] ; then
	update_owner=$3 # update owner
    else
	update_owner=0 # not update owner
    fi
    ##
    if [ $dir_mode == 0 ] ; then
	dmod='drwxrwx---'
	dmod_num=770
    else
	dmod='drwxr-x---'
	dmod_num=750
    fi
    ##
    st=(`stat -c "%A %U %G" "$infile"`) # drwxr-xr-x wangming yulab
    fmod=${st[0]}
    fown=${st[1]}
    fgrp=${st[2]}
    ##
    tag=0
    t1='pass'; t2='pass'; t3='pass'
    if [[ $fmod == d* ]] ; then
        # folder
	[[ $fmod != $dmod ]] && chmod $dmod_num $1 && t1=$fmod && tag=1
	[[ $fown != yuftp && $update_owner != 0 ]] && chown yuftp $1 && t2=$fown && tag=1
	[[ $fgrp != yulab ]] && chgrp yulab $1 && t3=$fgrp && tag=1
	[ $tag == 1 ] && echo -e "D: $t1\t$t2\t$t3\t$infile"
    elif [[ $fmod == -* ]] ; then
        # regular files
	[[ $fmod != -rw-r----- ]] && chmod 640 $1 && t1=$fmod && tag=1
	[[ $fown != yuftp && $update_owner != 0 ]] && chown yuftp $1 && t2=$fown && tag=1
	[[ $fgrp != yulab ]] && chgrp yulab $1 && t3=$fgrp && tag=1
	[ $tag == 1 ] && echo -e "F: $t1\t$t2\t$t3\t$infile"
    elif [[ $fmod == l* ]] ; then
        # symbolic link
        [[ $fmod != lrwxrwxrwx ]] && chmod 777 $1 && t1=$fmod && tag=1
	echo -e "L: $t1\t$t2\t$t3\t$infile"
    else
        # other type of files
	echo -e "X: File type not known" 
    fi
}

## backup files
function backup {
    # file
    # destination
    FROM=$1
    # dir1
    if [ x$2 != x ] ; then
	FROM_DIR=$2
    else
	FROM_DIR=$SHARE_DIR
    fi
    # dir2
    if [ x$3 != x ] ; then
	TO_DIR=$3
    else
	TO_DIR=$BACKUP_DIR
    fi
    # config
    dt=`date "+%Y%m%d_%H%M%S"`
    # check file existence in destination
    TO=${FROM/$FROM_DIR/$TO_DIR}
    if [ $FROM -nt $TO ] ; then
        ## FROM newer than CHECK, or CHECK not exists
        if [ -f $TO ] ; then
	    ## CHECK file older, rename by adding date to the tail
	    mv -f $TO ${TO}-${dt}.old
#	    echo $TO | rename -e 's/$/$dt/'
	fi
        ## copy FROM to TO (with directory structure)
        ## goto share_dir
        ## pwd
	PWD=$(readlink -f `pwd`)
	FROM_TMP=${FROM#$FROM_DIR\/}
	cd $FROM_DIR
        cp -rf -p --parents $FROM_TMP $TO_DIR
	cd $PWD
	echo "Backup: $FROM"
    else
	TAG=1
        #echo "Skip: $FROM"
    fi
}

## create log
DT_Y=`date "+%Y"`
DT_M=`date "+%m"`
DT_D=`date "+%d"`
DT_FULL=`date "+%Y%m%d"`
LOG_DIR="${BASE_DIR}/.log/${DT_Y}/${DT_M}"
LOG_FILE="${LOG_DIR}/${DT_FULL}.log"
LOG_TEMP="${LOG_DIR}/${DT_FULL}.temp"
LOG_ERR=${LOG_TEMP/.temp/.err}
[ ! -d $LOG_DIR ] && mkdir -p $LOG_DIR

###########
## start ##
###########
##1.update Share
## step1. search
find $SHARE_DIR > $LOG_TEMP 2> $LOG_ERR
## update files
while read LINE
do
    # replace blanks by "."
    echo "$LINE" | rename -e 's/ |\?/./g'
    NEW=`echo "$LINE" | sed -e 's/ |\?/./g'`
    # update files
    update $NEW 0 0 >> $LOG_FILE  # dir=770, owner=not change
done < $LOG_TEMP
## update Share_dir
chmod 770 $SHARE_DIR

## 2.backup files in SHARE_DIR
# single direction: Share -> Share_backup (hourly)
# > create file, not found in backup
# > update file (newer), in backup (rename old file by name)
find $SHARE_DIR > $LOG_TEMP 2> $LOG_ERR
while read LINE
do
    [[ $LINE == $SHARE_DIR ]] && continue
    backup $LINE >> $LOG_FILE
done < $LOG_TEMP
##

## 3. update backup
find $BACKUP_DIR > $LOG_TEMP 2> $LOG_ERR
while read LINE ; do
    update $LINE 1 1 >> $LOG_FILE # dir=750, owner=yuftp
done < $LOG_TEMP
##
rm $LOG_ERR $LOG_TEMP

##
MARK=`date "+%Y-%m-%d %H:%M:%S"`
echo "Finish @ $MARK"


###############################################################################
## Sync /data/yulab_ftp/Share/Papers/ to /data/yulab_ftp/For_intern/Papers/
###############################################################################
paperFrom="/data/yulab_ftp/Share/Papers"
paperTo="/data/yulab_ftp/For_intern/Papers"
rsync -aP ${paperFrom}/ ${paperTo}/ &>/dev/null
chown -R stu01:intern ${paperTo}

## EOF

### update folders that failed to access
#for i in $(seq 1 4) # go further 4-level in that directory
#do
#    DT_F=`date "+%Y-%m-%d %H:%M:%S"`
#    ERR_LINE=`wc -l $LOG_ERR | cut -d" " -f 1`
#    [ $ERR_LINE == 0 ] && echo "## Finish @ $DT_F" && break 
#    DIR_N=()
#    while read DL
#    do
#        DN=($DL)
#        DT=${DN%\'}
#        DT=${DT#\`}
#        DIR_N+=($DT)
#    done < $LOG_ERR
#    EXT_DD=$(IFS=" "; echo "${DIR_N[*]}")
#    find $EXT_DO > ${LOG_TEMP}.$i 2> $LOG_ERR
#    while read LINE
#    do
#        update $LINE
#    done < ${LOG_TEMP}.$i
#    rm ${LOG_TEMP}.$i $LOG_ERR
#done
#rm $LOG_TEMP

## EOF
