#!/bin/sh
# $1 SAM file
# $2 prefix for outputfile
~/BisPin/shell/split.sam.sh  $1 $2 0 74 &
~/BisPin/shell/split.sam.sh  $1 $2 75 125 &
~/BisPin/shell/split.sam.sh  $1 $2 126 1500 &
