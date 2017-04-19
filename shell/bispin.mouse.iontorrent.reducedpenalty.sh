#!/bin/sh

bash -c "time BisPin_align.py -x ~/BisPin/rescore_lower_gap_penalties.txt --gap_open=-200 --gap_extension=-20  -i 1 -I 2,3 -k -W /research/jsporter/Data/genome/GRCm38.p5/BisFAST/GRCm38.p5.fa /research/jsporter/Data/output/Mouse/BisFAST/Single/IonTorrent/SRR1534391.1mil.reduced_penalty.sam /research/jsporter/Data/reads/Mouse/IonTorrent/SRR1534391/SRR1534391.1mil.fastq" &> /research/jsporter/Data/output/Mouse/BisFAST/Single/IonTorrent/SRR1534391.1mil.reducedpenalty.out

bash -c "time BisPin_align.py -x ~/BisPin/rescore_lower_gap_penalties.txt --gap_open=-200 --gap_extension=-20  -i 1 -I 2,3 -k -W /research/jsporter/Data/genome/GRCm38.p5/BisFAST/GRCm38.p5.fa /research/jsporter/Data/output/Mouse/BisFAST/Single/IonTorrent/SRR1534392.1mil.reduced_penalty.sam /research/jsporter/Data/reads/Mouse/IonTorrent/SRR1534392/SRR1534392.1mil.fastq" &> /research/jsporter/Data/output/Mouse/BisFAST/Single/IonTorrent/SRR1534392.1mil.reducedpenalty.out
