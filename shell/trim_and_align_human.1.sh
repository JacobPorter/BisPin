#!/bin/sh
# $1 The reads directory and filename 
# $2 The output directory and filename
# $3 The number of threads

~/BisPin/shell/align_human.sh $1 $2.notrim $3

java -jar /home/jsporter/trimmomatic/classes/trimmomatic.jar SE -threads $3 -phred33 $1 $2.trimmomatic.fastq MAXINFO:25:0.5

~/BisPin/shell/align_human.sh $2.trimmomatic.fastq $2.trimmomatic.out $3

#User manual defaults
cutadapt -q 10 -o $2.cutadapt.fastq $1

~/BisPin/shell/align_human.sh $2.cutadapt.fastq $2.cutadapt.out $3

/home/jsporter/erne-2.1.1-linux/bin/erne-filter --query1 $1 --output-prefix $2.erne

~/BisPin/shell/align_human.sh $2.erne_1.fastq $2.erne.out $3 

#Dust is a measure of sequence complexity between 0 and 100.  A high score indicates a low complexity.
#~/Reaper/reaper -i $1 -basename $2.reaper -geom no-bc -3pa "" -tabu "" -dust-suffix 90 -qqq-check 50/15 --nozip 

#mv $2.reaper.lane.clean $2.reaper.fastq

#~/BisPin/shell/align_human.sh $2.reaper.fastq $2.reaper.out $3

#/home/jsporter/sickle/sickle se -f $1 -t sanger -o $2.sickle.fastq

#~/BisPin/shell/align_human.sh $2.sickle.fastq $2.sickle.out $3

#~/BisPin/InfoTrim.py -p $3  $1  1> $2.infotrim.fastq 2> $2.infotrim.info 

#~/BisPin/shell/align_human.sh $2.infotrim.fastq $2.infotrim.out $3
