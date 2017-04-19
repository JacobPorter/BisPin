#!/bin/sh
# $1 beginning iterate
# $2 ending iterate
# $3 number of reads to simulate
# $4 length of reads to simulate
for ((i=$1; i<=$2; i++)); do
   ~/BisPin/shell/simulation_paired.sh $i $3 $4
done
