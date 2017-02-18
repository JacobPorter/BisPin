#!/bin/sh

nohup bash -c "time BisPin_align.py -p 1  -k -i 1 -I 2 /research/jsporter/Data/genome/A_thaliana.TAIR10/BisPin/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa /research/jsporter/Data/output/A_thaliana/BisPin/Single/SRR5014638.500k.BisPin.sam  /research/jsporter/Data/reads/A_thaliana/Single/SRR5014638.500k.fastq" &> /research/jsporter/Data/output/A_thaliana/BisPin/Single/SRR5014638.500k.BisPin.out &

nohup bash -c "time BisPin_align.py -k -i 1 -I 2 /research/jsporter/Data/genome/A_thaliana.TAIR10/BisPin/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa /research/jsporter/Data/output/A_thaliana/BisPin/Paired-End/SRR4295457.500k.BisPin.sam  /research/jsporter/Data/reads/A_thaliana/Paired-End/SRR4295457_1.500k.fastq /research/jsporter/Data/reads/A_thaliana/Paired-End/SRR4295457_2.500k.fastq"  &> /research/jsporter/Data/output/A_thaliana/BisPin/Paired-End/SRR4295457.500k.BisPin.out &
