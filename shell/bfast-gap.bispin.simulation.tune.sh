#!/bin/sh

BisPin_align.py -n 7 -i 2 -k -W -N -g -b -10000000 -x /home/jsporter/bfast-gap-0.1.0/logistic_open_1.0_-15_constant_extension.txt /research/jsporter/Data/genome/GRCh38.p9/BisFAST/GRCh38.p9.fa /research/jsporter/Data/output/Human/BFAST-Gap/IonTorrent/Simulation/logistic/open_1.0_-15_extension_constant/reads.1m.dwgsim.human.ver0.1000024.logistic.bfast-gap.sam /research/jsporter/Data/reads/Human/IonTorrent/Simulation/reads.1m.dwgsim.human.ver0.1000024.fastq &> /research/jsporter/Data/output/Human/BFAST-Gap/IonTorrent/Simulation/logistic/open_1.0_-15_extension_constant/reads.1m.dwgsim.human.ver0.1000024.logistic.bfast-gap.out

cp /research/jsporter/Data/output/Human/BFAST-Gap/IonTorrent/Simulation/logistic/open_1.0_-15_extension_constant/*.a4.*.sam /research/jsporter/Data/output/Human/BFAST-Gap/IonTorrent/Simulation/logistic/open_1.0_-15_extension_constant/AUC/source/

gzip /research/jsporter/Data/output/Human/BFAST-Gap/IonTorrent/Simulation/logistic/open_1.0_-15_extension_constant/*.sam &

calculateAUCForBFAST_GAP.py -p 14 -o /research/jsporter/Data/output/Human/BFAST-Gap/IonTorrent/Simulation/logistic/open_1.0_-15_extension_constant/AUC/sink/   /research/jsporter/Data/genome/GRCh38.p9/BisFAST/GRCh38.p9.fa /research/jsporter/Data/output/Human/BFAST-Gap/IonTorrent/Simulation/logistic/open_1.0_-15_extension_constant/AUC/source/ /research/jsporter/Data/reads/Human/IonTorrent/Simulation/reads.1m.dwgsim.human.ver0.1000024.fastq 1000024 1> /research/jsporter/Data/output/Human/BFAST-Gap/IonTorrent/Simulation/logistic/open_1.0_-15_extension_constant/AUC/AUC.logistic.simulation.csv /research/jsporter/Data/output/Human/BFAST-Gap/IonTorrent/Simulation/logistic/open_1.0_-15_extension_constant/AUC/AUC.logistic.simulation.info
