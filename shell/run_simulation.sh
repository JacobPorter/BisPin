#!/bin/sh
# $1 beginning iterate
# $2 ending iterate
# $3 number of reads to simulate
for ((i=$1; i<=$2; i++)); do
   ~/BisPin/shell/simulation.sh $i $3
done
