#!/bin/sh
# 
nohup ~/DWGSIM/dwgsim -e 0.012 -E 0.012 -d 250 -s 30 -S 0 -N 1000000 -c 2 -1 200 -2 200 -f TACGTACGTCTGAGCATCGATCGATGTACAGC  /research/jsporter/Data/genome/GRCh38.p9/methyl-convert/GRCh38.p9.methyl-convert.GA.fa reads.1m.dwgsim.human.GA &> reads.1m.dwgsim.human.GA.out &

nohup ~/DWGSIM/dwgsim -e 0.012 -E 0.012 -d 250 -s 30 -S 0 -N 1000000 -c 2 -1 200 -2 200 -f TACGTACGTCTGAGCATCGATCGATGTACAGC  /research/jsporter/Data/genome/GRCh38.p9/methyl-convert/GRCh38.p9.methyl-convert.CT.fa reads.1m.dwgsim.human.CT &> reads.1m.dwgsim.human.CT.out &

