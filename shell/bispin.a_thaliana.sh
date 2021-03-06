#!/bin/sh

nohup bash -c "time BisFAST_align.py -p 1  -k -i 1 -I 2 /research/jsporter/Data/genome/A_thaliana.TAIR10/BisFAST/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa /research/jsporter/Data/output/A_thaliana/BisFAST/Single/SRR5014638.500k.bisfast.sam  /research/jsporter/Data/reads/A_thaliana/Single/SRR5014638.500k.fastq" &> /research/jsporter/Data/output/A_thaliana/BisFAST/Single/SRR5014638.500k.bisfast.out &

nohup bash -c "time BisFAST_align.py -k -i 1 -I 2 /research/jsporter/Data/genome/A_thaliana.TAIR10/BisFAST/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa /research/jsporter/Data/output/A_thaliana/BisFAST/Paired-End/SRR4295457.500k.bisfast.sam  /research/jsporter/Data/reads/A_thaliana/Paired-End/SRR4295457_1.500k.fastq /research/jsporter/Data/reads/A_thaliana/Paired-End/SRR4295457_2.500k.fastq"  &> /research/jsporter/Data/output/A_thaliana/BisFAST/Paired-End/SRR4295457.500k.bisfast.out &
