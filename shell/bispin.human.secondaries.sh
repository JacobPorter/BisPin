#!/bin/sh

BisPin_index.py -m 11110011110011110011 -w 14 -i 3 /research/jsporter/Data/genome/GRCh38.p9/BisFAST/GRCh38.p9.fa &> /research/jsporter/Data/genome/GRCh38.p9/BisFAST/3.index.out

BisPin_index.py -m 111100110110100111 -w 12 -i 4 /research/jsporter/Data/genome/GRCh38.p9/BisFAST/GRCh38.p9.fa &> /research/jsporter/Data/genome/GRCh38.p9/BisFAST/4.index.out

BisPin_index.py -m 111111111111 -w 12 -i 5 /research/jsporter/Data/genome/GRCh38.p9/BisFAST/GRCh38.p9.fa &> /research/jsporter/Data/genome/GRCh38.p9/BisFAST/5.index.out

BisPin_index.py -m 111111000111111 -w 12 -i 6 /research/jsporter/Data/genome/GRCh38.p9/BisFAST/GRCh38.p9.fa &> /research/jsporter/Data/genome/GRCh38.p9/BisFAST/6.index.out
