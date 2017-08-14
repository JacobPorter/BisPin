#!/bin/sh
nohup BisPin_align.py -n 6 -i 1 -I 2 -k -W -N -g -b -10000000 -x /home/jsporter/bfast-gap-0.1.0/logistic_scoring.txt /research/jsporter/Data/genome/GRCh38.p9/BisFAST/GRCh38.p9.fa /research/jsporter/Data/output/Human/BFAST-Gap/IonTorrent/Real/logistic/SRR2842546.1mil.bfast-gap.logistic.sam /research/jsporter/Data/reads/Human/IonTorrent/SRR2842546.1mil.fastq  &> /research/jsporter/Data/output/Human/BFAST-Gap/IonTorrent/Real/logistic/SRR2842546.1mil.bfast-gap.logistic.out 

nohup BisPin_align.py -n 6 -i 1 -I 2 -k -W -N -g -b -10000000 -x /home/jsporter/bfast-gap-0.1.0/logistic_scoring.txt /research/jsporter/Data/genome/GRCm38.p5/BisFAST/GRCm38.p5.fa  /research/jsporter/Data/output/Mouse/BFAST-Gap/IonTorrent/Real/logistic/SRR1534391.1mil.bfast-gap.logistic.sam /research/jsporter/Data/reads/Mouse/IonTorrent/SRR1534391/SRR1534391.1mil.fastq  &> /research/jsporter/Data/output/Mouse/BFAST-Gap/IonTorrent/Real/logistic/SRR1534391.1mil.bfast-gap.logistic.out
