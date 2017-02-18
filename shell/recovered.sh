#!/bin/sh
nohup BisPin_align.py -n 6  -i 1 -I 2,3 -k /research/jsporter/Data/genome/GRCm38.p5/BisFAST/GRCm38.p5.fa /research/jsporter/Data/output/Mouse/BisFAST/Hairpin/Single/E14-d0.L8-1.1mil.sam  /research/jsporter/Data/reads/Mouse/Hairpin/E14-d0.L8-1.1mil.fastq &> /research/jsporter/Data/output/Mouse/BisFAST/Hairpin/Single/E14-d0.L8-1.1mil.out &

nohup ./bismark --unmapped --ambiguous --sam --temp_dir /research/jsporter/Data/output/Mouse/bismark/Hairpin/Single/  -o /research/jsporter/Data/output/Mouse/bismark/Hairpin/Single/ /research/jsporter/Data/genome/GRCm38.p5/bismark/ /research/jsporter/Data/reads/Mouse/Hairpin/E14-d0.L8-1.1mil.fastq &> /research/jsporter/Data/output/Mouse/bismark/Hairpin/Single/E14-d0.L8.1mil.bismark.out &

nohup ./bowtie2 -x /research/jsporter/Data/genome/GRCm38.p5/bowtie2/GRCm38.p5.fa -U /research/jsporter/Data/reads/Mouse/Hairpin/BisFAST.hairpin.500k.R  -S /research/jsporter/Data/output/Mouse/bismark/Hairpin/Recovered/recovered.bowtie2.sam &> /research/jsporter/Data/output/Mouse/bismark/Hairpin/Recovered/recovered.bowtie2.out &

nohup ./bwameth.py --reference /research/jsporter/Data/genome/GRCm38.p5/bwa-meth/GRCm38.p5.fa /research/jsporter/Data/reads/Mouse/Hairpin/E14-d0.L8-1.1mil.fastq 1> /research/jsporter/Data/output/Mouse/bwa-meth/Hairpin/Single/E14-d0.L8-1.1mil.sam 2> /research/jsporter/Data/output/Mouse/bwa-meth/Hairpin/Single/E14-d0.L8-1.1mil.bwameth.out &

nohup ./bwa mem /research/jsporter/Data/genome/GRCm38.p5/bwa/GRCm38.p5.fa /research/jsporter/Data/reads/Mouse/Hairpin/BisFAST.hairpin.500k.R 1> /research/jsporter/Data/output/Mouse/bwa-meth/Hairpin/Recovered/recovered.bwa.sam 2> /research/jsporter/Data/output/Mouse/bwa-meth/Hairpin/Recovered/recovered.bwa.out &

nohup ./soap -u /research/jsporter/Data/output/Mouse/walt/Hairpin/Recovered/unmapped.soap.sam  -r 0  -D /research/jsporter/Data/genome/GRCm38.p5/soap/GRCm38.p5.fa.index -a /research/jsporter/Data/reads/Mouse/Hairpin/BisFAST.hairpin.500k.R.fq  -o /research/jsporter/Data/output/Mouse/walt/Hairpin/Recovered/recovered.soap.sam  &> /research/jsporter/Data/output/Mouse/walt/Hairpin/Recovered/recovered.soap.out &

nohup ./walt -i /research/jsporter/Data/genome/GRCm38.p5/walt/GRCm38.p5.multiLine.fa.index.walt.dbindex -r /research/jsporter/Data/reads/Mouse/Hairpin/E14-d0.L8-1.1mil.fastq -o /research/jsporter/Data/output/Mouse/walt/Hairpin/Single/E14-d0.L8-1.1mil.sam -t 4 -u -a &> /research/jsporter/Data/output/Mouse/walt/Hairpin/Single/walt.out &
