#!/bin/bash

###
# this script will backup the login history (from command `last`) 
#
# save record each day
DT=`date +%Y%m%d`
DY=`date +%d`
F="$HOME/bin/login_history/daybyday/${DT}.log"
last | awk -v DAY=$DY '{if($6==DAY) print $0}' > $F
#
# EOF
