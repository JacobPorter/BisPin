#!/bin/sh

nohup bash -c "time ~/bwa-meth/bwameth.py --reference /research/jsporter/Data/genome/A_thaliana.TAIR10/bwa-meth/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa /research/jsporter/Data/reads/A_thaliana/Paired-End/SRR4295457_1.500k.fastq /research/jsporter/Data/reads/A_thaliana/Paired-End/SRR4295457_2.500k.fastq" 1> /research/jsporter/Data/output/A_thaliana/bwa-meth/Paired-End/SRR4295457.500k.bwameth.sam 2> /research/jsporter/Data/output/A_thaliana/bwa-meth/Paired-End/SRR4295457.500k.bwameth.out &
