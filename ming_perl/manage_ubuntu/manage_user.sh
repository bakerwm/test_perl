#!/bin/bash

## This script describe how to manage users on Ubuntu 15.10 server


## Add user [id], create [HOME], assign to [group]
#sudo useradd -g [group] -m -s /bin/bash -d [HOME] [id]

## Assign password
#sudo passwd [id]

## Add user to exist [group]
#sudo usermod -G [group] [id]
#sudo gpasswd -a [id] [group]


## Lock user [id]
#sudo passwd -l [id]

## Unlock user [id]
#sudo passwd -u [id]


## Delete user and it's HOME
#sudo userdel --remove [id]

## Remove user from [group]
#sudo gpasswd -d [id] [group]


## Add FTP users (do not create home directory)
#sudo useradd -c "Real Name" -g yulab -G yulab -p passwd -s /bin/bash -M -d /data/yulab_ftp user

## Add SSH users (create /home/user)
#sudo useradd -c "Real Name" -g yulab -G yulab -p passwd -s /bin/bash -m -d /home/user user

##EOF
