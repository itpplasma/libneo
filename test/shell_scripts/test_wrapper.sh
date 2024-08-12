#!/bin/bash

# Short wrapper for the test functions.
# This is required, as an abort with 'stop' is still counted as sucess
# by ctest.
#
# Usage:
#   ./test_wrapper.sh test_to_execute
#
# Where 'test_to_execute' is the name of an executable.

test_to_execute=$1

if [ "x$2" == "xyes" ] ; then
  echo "running test_arnoldi.py"
  pwd
  ../test/python_scripts/test_arnoldi.py
fi

./$test_to_execute | grep -i " STOP "
if [ "x$?" == "x0" ] ; then
  exit 1
else
  exit 0
fi
