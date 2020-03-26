#!/bin/bash
SCRIPTNAME=$(basename $0 .sh)
EXIT_SUCCESS=0
EXIT_FAILURE=1

function usage {
 echo "Usage: ./$SCRIPTNAME [-d (dry-run)]" >&2
 [[ $# -eq 1 ]] && exit $1 || exit $EXIT_FAILURE
}

DRY=0
while getopts 'dh' OPTION ; do
 case $OPTION in
 d) echo "Running in dry-mode, no job is started"; DRY=1
 ;;
 h) usage $EXIT_SUCCESS
 ;;
 \?) echo "Unknown Option \"-$OPTARG\"." >&2
 usage $EXIT_ERROR
 ;;
 esac
done

for i in $(cat jobs_list.txt | grep -v '^#' | grep -v "^$"); do
  echo "Queueing $i..."
  if [ $DRY -eq 0 ] ; then
  	cd $i
  	condor_submit condor.submit
  	cd ..
  fi
done;
