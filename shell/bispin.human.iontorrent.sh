#!/bin/sh
# $1 The reads filename.
# $2 The number of threads for each process

bash -c "time BisPin_align.py -n $2 -i 1 -I 2 -k -W /research/jsporter/Data/genome/GRCh38.p9/BisFAST/GRCh38.p9.fa /research/jsporter/Data/output/Human/BisFAST/Single/IonTorrent/$1.sam /research/jsporter/Data/reads/Human/IonTorrent/$1" &> /research/jsporter/Data/output/Human/BisFAST/Single/IonTorrent/$1.bispin.out
