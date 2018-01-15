#!/bin/sh

BisPin_align.py -n 6 -i 2 -k -W -N -g -b -10000000 -x /home/jsporter/bfast-gap-0.1.0/logistic_open_1.0_-15_constant_extension.txt /research/jsporter/Data/genome/GRCh38.p9/BisFAST/GRCh38.p9.fa /research/jsporter/Data/output/Human/BFAST-Gap/IonTorrent/Real/logistic/open_1.0_-15_extension_constant/SRR2842546.1mil.bispin.logistic.sam /research/jsporter/Data/reads/Human/IonTorrent/SRR2842546.1mil.fastq &> /research/jsporter/Data/output/Human/BFAST-Gap/IonTorrent/Real/logistic/open_1.0_-15_extension_constant/SRR2842546.1mil.bfast-gap.out

cp /research/jsporter/Data/output/Human/BFAST-Gap/IonTorrent/Real/logistic/open_1.0_-15_extension_constant/*.a4.*.sam /research/jsporter/Data/output/Human/BFAST-Gap/IonTorrent/Real/logistic/open_1.0_-15_extension_constant/AUC/source/

gzip /research/jsporter/Data/output/Human/BFAST-Gap/IonTorrent/Real/logistic/open_1.0_-15_extension_constant/*.sam &

calculateAUCForBFAST_GAP.py -p 13 -o /research/jsporter/Data/output/Human/BFAST-Gap/IonTorrent/Real/logistic/open_1.0_-15_extension_constant/AUC/sink/   /research/jsporter/Data/genome/GRCh38.p9/BisFAST/GRCh38.p9.fa /research/jsporter/Data/output/Human/BFAST-Gap/IonTorrent/Real/logistic/open_1.0_-15_extension_constant/AUC/source/ /research/jsporter/Data/reads/Human/IonTorrent/Real/SRR2842546.1mil.fastq 1000000 1> /research/jsporter/Data/output/Human/BFAST-Gap/IonTorrent/Real/logistic/open_1.0_-15_extension_constant/AUC/AUC.logistic.csv /research/jsporter/Data/output/Human/BFAST-Gap/IonTorrent/Real/logistic/open_1.0_-15_extension_constant/AUC/AUC.logistic.info
