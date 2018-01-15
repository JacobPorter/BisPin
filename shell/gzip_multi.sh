#!/bin/sh
# Launches multiple gzip processing cores.
# $1 file to gzip
# $2 number of cores to use
# CORES=$(grep -c '^processor' /proc/cpuinfo)
CORES=$2
find $1  -type f -print0 | xargs -0 -n 1 -P $CORES gzip
