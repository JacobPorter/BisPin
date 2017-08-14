#!/bin/sh
# $1 reads file
# $2 number of threads
BisPin_align.py -n $2 -i 1 -I 2,3,4,5,6 -W /research/jsporter/Data/genome/GRCh38.p9/BisFAST/GRCh38.p9.fa /research/jsporter/Data/output/Human/BisFAST/Single/IonTorrent/$1.bispin.sam $1 &> /research/jsporter/Data/output/Human/BisFAST/Single/IonTorrent/$1.bispin.multiIndexes.out 
