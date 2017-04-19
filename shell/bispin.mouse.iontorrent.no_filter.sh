#!/bin/sh

bash -c "time BisPin_align.py -n 4  -m -1000  -i 1 -I 2,3 -k -W /research/jsporter/Data/genome/GRCm38.p5/BisFAST/GRCm38.p5.fa /research/jsporter/Data/output/Mouse/BisFAST/Single/IonTorrent/SRR1534391.1mil.no_filter.sam /research/jsporter/Data/reads/Mouse/IonTorrent/SRR1534391/SRR1534391.1mil.fastq" &> /research/jsporter/Data/output/Mouse/BisFAST/Single/IonTorrent/SRR1534391.1mil.no_filter.out

bash -c "time BisPin_align.py -n 4  -m -1000  -i 1 -I 2,3 -k -W /research/jsporter/Data/genome/GRCm38.p5/BisFAST/GRCm38.p5.fa /research/jsporter/Data/output/Mouse/BisFAST/Single/IonTorrent/SRR1534392.1mil.no_filter.sam /research/jsporter/Data/reads/Mouse/IonTorrent/SRR1534392/SRR1534392.1mil.fastq" &> /research/jsporter/Data/output/Mouse/BisFAST/Single/IonTorrent/SRR1534392.1mil.no_filter.out

