#!/bin/sh
# $1 read min length
# $2 read max length
# $3 read file name
# $4 number of reads
# $5 threads to use
NUMREADS=`echo "$4 * 4" | bc`
NEWFILE=$3.$1-$2
NEWNUMFILE=$3.$1-$2.$4
lengthSelector.py -t FASTQ $1 $2 $3 1> $NEWFILE.fastq 2> $NEWFILE.select.out
head -n $NUMREADS $NEWFILE.fastq 1> $NEWNUMFILE.fastq 2> $NEWNUMFILE.out
rm $NEWFILE.fastq

~/bismark_v0.16.3/bismark -p $5  --ambiguous --unmapped  --sam --temp_dir /research/jsporter/Data/output/Human/bismark/Single/IonTorrent/ -o /research/jsporter/Data/output/Human/bismark/Single/IonTorrent/$NEWNUMFILE/  /research/jsporter/Data/genome/GRCh38.p9/bismark/ $NEWNUMFILE.fastq  &> /research/jsporter/Data/output/Human/bismark/Single/IonTorrent/$NEWNUMFILE.bismark.out

~/bwa-meth/bwameth.py --reference /research/jsporter/Data/genome/GRCh38.p9/bwa-meth/GRCh38.p9.fa $NEWNUMFILE.fastq 1> /research/jsporter/Data/output/Human/bwa-meth/Single/IonTorrent/$NEWNUMFILE.bwameth.sam 2> /research/jsporter/Data/output/Human/bwa-meth/Single/IonTorrent/$NEWNUMFILE.bwameth.out

countBWA.py /research/jsporter/Data/output/Human/bwa-meth/Single/IonTorrent/$NEWNUMFILE.bwameth.sam 1> /research/jsporter/Data/output/Human/bwa-meth/Single/IonTorrent/$NEWNUMFILE.bwameth.count 2> /research/jsporter/Data/output/Human/bwa-meth/Single/IonTorrent/$NEWNUMFILE.bwameth.count.out

~/walt-1.0/bin/walt -t $5 -i /research/jsporter/Data/genome/GRCh38.p9/walt/GRCh38.p9.multiLine.fa.index.walt.dbindex -o /research/jsporter/Data/output/Human/walt/Single/IonTorrent/$NEWNUMFILE.walt.sam -r $NEWNUMFILE.fastq  &> /research/jsporter/Data/output/Human/walt/Single/IonTorrent/$NEWNUMFILE.walt.out 

/home/jsporter/tabsat/tools/bismark_tmap/bismark --path_to_program ~/TMAP/ --temp_dir /research/jsporter/Data/output/Human/ --tmap -o /research/jsporter/Data/output/Human/tabsat/IonTorrent/$NEWNUMFILE/ /research/jsporter/Data/genome/GRCh38.p9/tabsat/ $NEWNUMFILE.fastq &> /research/jsporter/Data/output/Human/tabsat/IonTorrent/$NEWNUMFILE.tabsat.out

BisPin_align.py -n $5 -i 1 -I 2 -W /research/jsporter/Data/genome/GRCh38.p9/BisFAST/GRCh38.p9.fa /research/jsporter/Data/output/Human/BisFAST/Single/IonTorrent/$NEWNUMFILE.bispin.sam $NEWNUMFILE.fastq &> /research/jsporter/Data/output/Human/BisFAST/Single/IonTorrent/$NEWNUMFILE.bispin.out

alignmentScoreDistribution.py -b 20 /research/jsporter/Data/output/Human/BisFAST/Single/IonTorrent/$NEWNUMFILE.bispin.sam 1> /research/jsporter/Data/output/Human/BisFAST/Single/IonTorrent/$NEWNUMFILE.bispin.score 2> /research/jsporter/Data/output/Human/BisFAST/Single/IonTorrent/$NEWNUMFILE.bispin.score.out 

rm $NEWNUMFILE.fastq
