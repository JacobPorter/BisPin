#!/bin/sh

nohup bash -c "time ~/walt-1.0/bin/walt -a -u  -i /research/jsporter/Data/genome/A_thaliana.TAIR10/walt/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.walt.dbindex  -o /research/jsporter/Data/output/A_thaliana/walt/Paired-End/SRR4295457.500k.walt.sam -1 /research/jsporter/Data/reads/A_thaliana/Paired-End/SRR4295457_1.500k.fastq  -2 /research/jsporter/Data/reads/A_thaliana/Paired-End/SRR4295457_2.500k.fastq" &> /research/jsporter/Data/output/A_thaliana/walt/Paired-End/SRR4295457.500k.walt.out &
