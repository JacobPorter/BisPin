#!/bin/sh
nohup bash -c "time BisPin_align.py -p 1  -i 1 -I 2  -k /research/jsporter/Data/genome/GRCh38.p9/BisPin/GRCh38.p9.fa /research/jsporter/Data/output/Human/BisPin/Single/Nondirectional/SRR342552.500k.BisPin.sam /research/jsporter/Data/reads/Human/Single/Nondirectional/SRR342552.500k.fastq" &> /research/jsporter/Data/output/Human/BisPin/Single/Nondirectional/SRR342552.500k.BisPin.out &

nohup bash -c "time BisPin_align.py -p 2  -i 1 -I 2  -k /research/jsporter/Data/genome/GRCh38.p9/BisPin/GRCh38.p9.fa /research/jsporter/Data/output/Human/BisPin/Paired-End/PBAT/SRR2013793.500k.BisPin.sam /research/jsporter/Data/reads/Human/Paired-End/PBAT/SRR2013793_1.500k.fastq /research/jsporter/Data/reads/Human/Paired-End/PBAT/SRR2013793_2.500k.fastq" &> /research/jsporter/Data/output/Human/BisPin/Paired-End/PBAT/SRR2013793.500k.BisPin.out &
