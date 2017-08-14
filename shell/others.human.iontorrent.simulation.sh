#!/bin/sh
# $1 number of threads
# $2 reads filename
# $3 the number of reads
~/bismark_v0.16.3/bismark -p $1  --ambiguous --unmapped  --sam --temp_dir /research/jsporter/Data/output/Human/bismark/Single/IonTorrent/Simulation/ -o /research/jsporter/Data/output/Human/bismark/Single/IonTorrent/Simulation/$2/  /research/jsporter/Data/genome/GRCh38.p9/bismark/ /research/jsporter/Data/reads/Human/IonTorrent/Simulation/$2  &> /research/jsporter/Data/output/Human/bismark/Single/IonTorrent/Simulation/$2.bismark.out

~/bwa-meth/bwameth.py --reference /research/jsporter/Data/genome/GRCh38.p9/bwa-meth/GRCh38.p9.fa /research/jsporter/Data/reads/Human/IonTorrent/Simulation/$2 1> /research/jsporter/Data/output/Human/bwa-meth/Single/IonTorrent/Simulation/$2.bwameth.sam 2> /research/jsporter/Data/output/Human/bwa-meth/Single/IonTorrent/Simulation/$2.bwameth.err

countBWA.py /research/jsporter/Data/output/Human/bwa-meth/Single/IonTorrent/Simulation/$2.bwameth.sam 1> /research/jsporter/Data/output/Human/bwa-meth/Single/IonTorrent/Simulation/$2.bwameth.count 2> /research/jsporter/Data/output/Human/bwa-meth/Single/IonTorrent/Simulation/$2.bwameth.count.err

calculateSimulationAccuracy.py -d /research/jsporter/Data/output/Human/bwa-meth/Single/IonTorrent/Simulation/$2.bwameth.sam $3 1> /research/jsporter/Data/output/Human/bwa-meth/Single/IonTorrent/Simulation/$2.bwameth.acc 2> /research/jsporter/Data/output/Human/bwa-meth/Single/IonTorrent/Simulation/$2.bwameth.acc.err

~/walt-1.0/bin/walt -t $1 -i /research/jsporter/Data/genome/GRCh38.p9/walt/GRCh38.p9.multiLine.fa.index.walt.dbindex -o /research/jsporter/Data/output/Human/walt/Single/IonTorrent/Simulation/$2.walt.sam -r /research/jsporter/Data/reads/Human/IonTorrent/Simulation/$2  &> /research/jsporter/Data/output/Human/walt/Single/IonTorrent/Simulation/$2.walt.out

calculateSimulationAccuracy.py -d /research/jsporter/Data/output/Human/walt/Single/IonTorrent/Simulation/$2.walt.sam $3 1> /research/jsporter/Data/output/Human/walt/Single/IonTorrent/Simulation/$2.walt.acc 2> /research/jsporter/Data/output/Human/walt/Single/IonTorrent/Simulation/$2.walt.acc.err
