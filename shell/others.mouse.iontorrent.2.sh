#!/bin/sh
~/bismark_v0.16.3/bismark -p 2  --ambiguous --unmapped  --sam --temp_dir /research/jsporter/Data/output/Mouse/bismark/Single/IonTorrent/ -o /research/jsporter/Data/output/Mouse/bismark/Single/IonTorrent/SRR1534392/  /research/jsporter/Data/genome/GRCm38.p5/bismark/ /research/jsporter/Data/reads/Mouse/IonTorrent/SRR1534392/SRR1534392.1mil.fastq  &> /research/jsporter/Data/output/Mouse/bismark/Single/IonTorrent/SRR1534392.1mil.bismark.out

~/bwa-meth/bwameth.py --reference /research/jsporter/Data/genome/GRCm38.p5/bwa-meth/GRCm38.p5.fa /research/jsporter/Data/reads/Mouse/IonTorrent/SRR1534392/SRR1534392.1mil.fastq 1> /research/jsporter/Data/output/Mouse/bwa-meth/Single/IonTorrent/SRR1534392.1mil.bwameth.sam 2> /research/jsporter/Data/output/Mouse/bwa-meth/Single/IonTorrent/SRR1534392.1mil.bwameth.err

~/walt-1.0/bin/walt -t 2 -i /research/jsporter/Data/genome/GRCm38.p5/walt/GRCm38.p5.multiLine.fa.index.walt.dbindex -o /research/jsporter/Data/output/Mouse/walt/Single/IonTorrent/SRR1534392.1mil.walt.sam -r /research/jsporter/Data/reads/Mouse/IonTorrent/SRR1534392/SRR1534392.1mil.fastq  &> /research/jsporter/Data/output/Mouse/walt/Single/IonTorrent/SRR1534392.1mil.walt.out

