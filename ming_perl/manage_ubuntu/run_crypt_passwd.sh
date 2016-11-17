#!/bin/bash

while read LINE 
do
    a=($LINE)
    CRYPT=`openssl passwd -1 -salt xyz ${a[2]}`
    echo "${a[0]} ${a[1]} $CRYPT"
done < user_id_list.txt
