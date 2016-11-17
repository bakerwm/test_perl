#!/bin/bash

# Create users-specific directories in /data/yulab_ftp/01.Members/
#

## 1. check current users
## READ current users
cat /etc/passwd | cut -d: -f 1  > /home/wangming/bin/user_id_current.txt

## Create new FTP dirs for new users
while IFS= read -r LINE
do
    ## FTP directory
    USER_FTP="/data/yulab_ftp/01.Members/$LINE"
    if [[ ! -d $USER_FTP ]] ; then
        echo "Create FTP dir for: $LINE [ $USER_FTP ]" 
        sudo mkdir -p $USER_FTP
        ## change the owership
        USER_uid=`id -u $LINE`
        USER_gid=`id -g $LINE`
        sudo chown $USER_uid:$USER_gid $USER_FTP
    fi
done < <(cat /home/wangming/bin/user_id_fixed.txt /home/wangming/bin/user_id_current.txt | sort | uniq -u)

## EOF
