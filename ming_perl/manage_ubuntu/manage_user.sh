#!/bin/bash

## This script describe how to manage users on Ubuntu 15.10 server


## Add user [id], create [HOME], assign to [group]
#sudo useradd -g [group] -m -s /bin/bash -d [HOME] [user]

## Assign password
#sudo passwd [user]

## Add user to exist [group]
#sudo usermod -G [group] [user]
#sudo gpasswd -a [user] [group]


## Lock user [user]
#sudo passwd -l [user]

## Unlock user [user]
#sudo passwd -u [user]


## Delete user and it's HOME
#sudo userdel --remove [user]

## Remove user from [group]
#sudo gpasswd -d [user] [group]


## Add FTP users (do not create home directory)
#sudo useradd -c "Real Name" -g yulab -G yulab -p passwd -s /bin/bash -M -d /data/yulab_ftp user

## Add SSH users (create /home/user)
#sudo useradd -c "Real Name" -g yulab -G yulab -p passwd -s /bin/bash -m -d /home/user user

##EOF
