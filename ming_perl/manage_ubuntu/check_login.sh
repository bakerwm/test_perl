#!/bin/bash

###
# this script will backup the login history (from command `last`) 
#
# save record each day
outdir="${HOME}/bin/login_history/byday"
dt=$(date +%Y%m%d)
year=$(date +%Y)
mon=$(date +%m)
day=$(date +%d)

##
subdir="${outdir}/${year}/${mon}"
log="${subdir}/${dt}.last.cmd.log"
mkdir -p ${subdir}
last | awk -v d=${day} '{if($6==d) print $0}' > ${log}

#DT=`date +%Y%m%d`
#DY=`date +%d`
##DT="2016120$1"
##DY=$1
#F="$HOME/bin/login_history/daybyday/${DT}.log"
#last | awk -v DAY=$DY '{if($6==DAY) print $0}' #> $F
##
## EOF
