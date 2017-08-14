#!/bin/sh
# $1 number of processes
# $2 input reads file name
# $3 total reads in the file
nohup BisPin_align.py -p 0 -n $1  -i 1 -I 2 -k /research/jsporter/Data/genome/GRCh38.p9/BisFAST/GRCh38.p9.fa /research/jsporter/Data/output/Human/BisFAST/Single/IonTorrent/Simulation/$2.bispin.default.sam /research/jsporter/Data/reads/Human/IonTorrent/Simulation/$2 &> /research/jsporter/Data/output/Human/BisFAST/Single/IonTorrent/Simulation/$2.bispin.default.out

nohup calculateSimulationAccuracy.py -d /research/jsporter/Data/output/Human/BisFAST/Single/IonTorrent/Simulation/$2.bispin.default.sam $3 1> /research/jsporter/Data/output/Human/BisFAST/Single/IonTorrent/Simulation/$2.bispin.default.acc 2> /research/jsporter/Data/output/Human/BisFAST/Single/IonTorrent/Simulation/$2.bispin.default.acc.err

nohup BisPin_align.py -p 0 -n $1 -m -20000 -i 1 -I 2 -k /research/jsporter/Data/genome/GRCh38.p9/BisFAST/GRCh38.p9.fa /research/jsporter/Data/output/Human/BisFAST/Single/IonTorrent/Simulation/$2.bispin.nofilter.sam /research/jsporter/Data/reads/Human/IonTorrent/Simulation/$2 &> /research/jsporter/Data/output/Human/BisFAST/Single/IonTorrent/Simulation/$2.bispin.nofilter.out

nohup calculateSimulationAccuracy.py -d /research/jsporter/Data/output/Human/BisFAST/Single/IonTorrent/Simulation/$2.bispin.nofilter.sam $3 1> /research/jsporter/Data/output/Human/BisFAST/Single/IonTorrent/Simulation/$2.bispin.nofilter.acc 2> /research/jsporter/Data/output/Human/BisFAST/Single/IonTorrent/Simulation/$2.bispin.nofilter.acc.err

nohup BisPin_align.py -p 0 -n $1 -m -20000 --gap_open=-200 --gap_extension=-20 -i 1 -I 2 /research/jsporter/Data/genome/GRCh38.p9/BisFAST/GRCh38.p9.fa /research/jsporter/Data/output/Human/BisFAST/Single/IonTorrent/Simulation/$2.bispin.nofilter_reducedpenalty.sam /research/jsporter/Data/reads/Human/IonTorrent/Simulation/$2 &> /research/jsporter/Data/output/Human/BisFAST/Single/IonTorrent/Simulation/$2.bispin.nofilter_reducedpenalty.out

nohup calculateSimulationAccuracy.py -d /research/jsporter/Data/output/Human/BisFAST/Single/IonTorrent/Simulation/$2.bispin.nofilter_reducedpenalty.sam $3 1> /research/jsporter/Data/output/Human/BisFAST/Single/IonTorrent/Simulation/$2.bispin.nofilter_reducedpenalty.acc 2> /research/jsporter/Data/output/Human/BisFAST/Single/IonTorrent/Simulation/$2.bispin.nofilter_reducedpenalty.acc.err
