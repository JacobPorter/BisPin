#!/bin/sh
# $1 number of reads
# $2 number of processes
# $3 read file name
~/BisPin/shell/lengthselector.human.sh 1 75 $3 $1 $2 
~/BisPin/shell/lengthselector.human.sh 76 150 $3 $1 $2 
~/BisPin/shell/lengthselector.human.sh 151 225 $3 $1 $2 
~/BisPin/shell/lengthselector.human.sh 226 300 $3 $1 $2 
~/BisPin/shell/lengthselector.human.sh 301 375 $3 $1 $2 
