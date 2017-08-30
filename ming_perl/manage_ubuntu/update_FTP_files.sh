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
# 2017-06-14 version 0.4
# Illustrate the workflow
# 1. backup directories:
# from: /data/yulab_ftp/Share
# to:   /data/yulab_ftp/Archive/Share_backup
# 
# 2. change permissions
#   a) Share/
#      dirs: 770
#      regular files: 640
#      symbolic links: 777
#   b) Share.backup/
#      ower: yuftp:yulab
#      dirs: 750
#      regular files: 640
#      symbolic links: 777
#   c) other folders/
#      ower: yuftp:yulab
#      dirs: 750
#      regular files: 640
#      symbolic links: 777
# 2016-12-02 version 0.3
#   - add backup function
#   - hourly update (*:59)
# 2016-10-26 version 0.2
#   - rewrite update function
#   - daily update (4:00 am)

##===============##
## global config ##
##===============##
ftpdir="/data/yulab_ftp"
share="${ftpdir}/Share"
shareBackup="${ftpdir}/Archive/Share_backup"
logdir="${HOME}/work/bin/temp/manage_ubuntu/logs"

## 1. rename files with blanks in name 
function renameFiles {
    # replace the blanks in the file name
    indir=$1
    [[ ! -d ${indir} ]] && echo "dir not exists - [${indir}]" && exit 1
    # list all files
    tmpfile="/tmp/$$.update_files.list"
    find ${indir} -type f > ${tmpfile}
    while read line ; do 
        dir=$(dirname "${line}")
        prefix=$(basename "${line}")
        prefixnew=${prefix// /.}
        fnew="${dir}/${prefixnew}"
        if [[ "${line}" != ${fnew} ]] ; then
            mv -f "${line}" ${fnew}
            echo "rename: [${line}] [${fnew}]"
        fi
    done < ${tmpfile}
    rm ${tmpfile}
}

function renameDirs {
    # replace the blanks in the dirname
    indir=$1
    [[ ! -d ${indir} ]] && echo "dir not exists - [${indir}]" && exit 1
    # move files to new dir in each depth
    for n in {1..10} ; do # suppose max 10 levels of depth
        # list all dirs
        tmpdir="/tmp/$$.update_dirs.list"
        find ${indir} -maxdepth ${n} -type d > ${tmpdir}
        while read line ; do 
            dir=$(dirname "${line}")
            prefix=$(basename "${line}")
            prefixnew=${prefix// /.}
            dnew="${dir}/${prefixnew}"
            if [[ "${line}" != ${dnew} ]] ; then
                mv -f "${line}" ${dnew}
                echo "rename: [${line}] [${dnew}]"
            fi
        done < ${tmpdir}
        rm ${tmpdir}
    done
}

## 1. rename files in Share/
renameDirs ${share}
renameFiles ${share}

## 2. update owership of files in Share/
find ${share} -type d | xargs chmod 777
find ${share} -type f | xargs chmod 640 
[[ $(find ${share} -type l | wc -l) -gt 0 ]] && find ${share} -type l | xargs chmod 777

## 3. backup Share/
rsync -avq --exclude "Old_files" --max-size=100m ${share}/ ${shareBackup}/ 

## 4. update owership and permission in Share_backup/
find ${shareBackup} -type d | xargs chmod 750
find ${shareBackup} -type f | xargs chmod 640
[[ $(find ${shareBackup} -type l | wc -l) -gt 0 ]] && find ${shareBackup} -type l | xargs chmod 777
find ${shareBackup} | xargs chown yuftp:yulab 

################################################################################
### Sync /data/yulab_ftp/Share/Papers/ to /data/yulab_ftp/For_intern/Papers/
################################################################################
rsync -avq /data/yulab_ftp/Share/Papers/ /data/yulab_ftp/For_intern/Papers/
chown -R stu01:intern /data/yulab_ftp/For_intern/Papers/

echo "Finish @ $(date "+%Y-%m-%d %H:%M:%S")"

## EOF
