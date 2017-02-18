#!/bin/sh
nohup bash -c "time ~/bismark_v0.16.3/bismark --non-directional --ambiguous --unmapped --sam --temp_dir /research/jsporter/Data/output/A_thaliana/bismark/Single/ -o /research/jsporter/Data/output/A_thaliana/bismark/Single/ /research/jsporter/Data/genome/A_thaliana.TAIR10/bismark/ /research/jsporter/Data/reads/A_thaliana/Single/SRR5014638.500k.fastq" &>  /research/jsporter/Data/output/A_thaliana/bismark/Single/SRR5014638.500k.bismark.out &

nohup bash -c "time ~/bismark_v0.16.3/bismark --ambiguous --unmapped --sam --temp_dir /research/jsporter/Data/output/A_thaliana/bismark/Paired-End/ -o /research/jsporter/Data/output/A_thaliana/bismark/Paired-End/ /research/jsporter/Data/genome/A_thaliana.TAIR10/bismark/ -1 /research/jsporter/Data/reads/A_thaliana/Paired-End/SRR4295457_1.500k.fastq -2 /research/jsporter/Data/reads/A_thaliana/Paired-End/SRR4295457_2.500k.fastq" &> /research/jsporter/Data/output/A_thaliana/bismark/Paired-End/SRR4295457.500k.bismark.out &


