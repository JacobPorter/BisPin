#!/bin/sh
# $1 The reads filename

~/bismark_v0.16.3/bismark -p 2  --ambiguous --unmapped  --sam --temp_dir /research/jsporter/Data/output/Human/bismark/Single/IonTorrent/ -o /research/jsporter/Data/output/Human/bismark/Single/IonTorrent/$1/  /research/jsporter/Data/genome/GRCh38.p9/bismark/ /research/jsporter/Data/reads/Human/IonTorrent/$1  &> /research/jsporter/Data/output/Human/bismark/Single/IonTorrent/$1.bismark.out

~/bwa-meth/bwameth.py --reference /research/jsporter/Data/genome/GRCh38.p9/bwa-meth/GRCh38.p9.fa /research/jsporter/Data/reads/Human/IonTorrent/$1 1> /research/jsporter/Data/output/Human/bwa-meth/Single/IonTorrent/$1.bwameth.sam 2> /research/jsporter/Data/output/Human/bwa-meth/Single/IonTorrent/$1.bwameth.err

~/walt-1.0/bin/walt -t 2 -i /research/jsporter/Data/genome/GRCh38.p9/walt/GRCh38.p9.multiLine.fa.index.walt.dbindex -o /research/jsporter/Data/output/Human/walt/Single/IonTorrent/$1.walt.sam -r /research/jsporter/Data/reads/Human/IonTorrent/$1  &> /research/jsporter/Data/output/Human/walt/Single/IonTorrent/$1.walt.out

