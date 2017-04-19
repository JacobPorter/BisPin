#!/bin/sh
# $1 read iterate number
# $2 number of reads to simulate
~/Sherman_v0.1.7/Sherman -n $2  -e 2 -l 150 -CG 20 -CH 2 --genome_folder /research/jsporter/Data/genome/GRCm38.p5/genome_folder/ &> /research/jsporter/Data/reads/Mouse/Single/Simulation/simulated.single.$1.out

mv simulated.fastq /research/jsporter/Data/reads/Mouse/Single/Simulation/simulated.single.$1.fastq

BisPin_align.py -W -i 1 -I 2,3 -n 6 -u /research/jsporter/Data/genome/GRCm38.p5/BisFAST/GRCm38.p5.fa /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.single.$1.sam /research/jsporter/Data/reads/Mouse/Single/Simulation/simulated.single.$1.fastq &> /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.single.$1.bispin.out

calculateSimulationAccuracy.py /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.single.$1.sam $2 1> /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.single.$1.acc 2> /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.single.$1.acc.out

~/bismark_v0.16.3/bismark -p 6  --ambiguous --unmapped  --sam --temp_dir /research/jsporter/Data/output/Mouse/bismark/Single/Simulation/ -o /research/jsporter/Data/output/Mouse/bismark/Single/Simulation/$1/  /research/jsporter/Data/genome/GRCm38.p5/bismark/ /research/jsporter/Data/reads/Mouse/Single/Simulation/simulated.single.$1.fastq &> /research/jsporter/Data/output/Mouse/bismark/Single/Simulation/simulated.single.$1.bismark.out

calculateSimulationAccuracy.py /research/jsporter/Data/output/Mouse/bismark/Single/Simulation/$1/simulated.single.$1_bismark_bt2.sam $2 1> /research/jsporter/Data/output/Mouse/bismark/Single/Simulation/$1/simulated.single.$1.acc 2> /research/jsporter/Data/output/Mouse/bismark/Single/Simulation/$1/simulated.single.$1.acc.out

mv /research/jsporter/Data/output/Mouse/bismark/Single/Simulation/$1/simulated.single.$1.acc /research/jsporter/Data/output/Mouse/bismark/Single/Simulation/simulated.single.$1.acc

~/bwa-meth/bwameth.py --reference /research/jsporter/Data/genome/GRCm38.p5/bwa-meth/GRCm38.p5.fa /research/jsporter/Data/reads/Mouse/Single/Simulation/simulated.single.$1.fastq 1> /research/jsporter/Data/output/Mouse/bwa-meth/Single/Simulation/simulated.single.$1.sam 2> /research/jsporter/Data/output/Mouse/bwa-meth/Single/Simulation/simulated.single.$1.bwameth.err

calculateSimulationAccuracy.py /research/jsporter/Data/output/Mouse/bwa-meth/Single/Simulation/simulated.single.$1.sam $2 1> /research/jsporter/Data/output/Mouse/bwa-meth/Single/Simulation/simulated.single.$1.acc 2> /research/jsporter/Data/output/Mouse/bwa-meth/Single/Simulation/simulated.single.$1.acc.out

~/walt-1.0/bin/walt -t 6 -i /research/jsporter/Data/genome/GRCm38.p5/walt/GRCm38.p5.multiLine.fa.index.walt.dbindex -o /research/jsporter/Data/output/Mouse/walt/Single/Simulation/simulated.single.$1.sam -r /research/jsporter/Data/reads/Mouse/Single/Simulation/simulated.single.$1.fastq &> /research/jsporter/Data/output/Mouse/walt/Single/Simulation/simulated.single.$1.walt.out

calculateSimulationAccuracy.py /research/jsporter/Data/output/Mouse/walt/Single/Simulation/simulated.single.$1.sam $2 1> /research/jsporter/Data/output/Mouse/walt/Single/Simulation/simulated.single.$1.acc 2> /research/jsporter/Data/output/Mouse/walt/Single/Simulation/simulated.single.$1.acc.out

rm /research/jsporter/Data/reads/Mouse/Single/Simulation/simulated.single.$1.fastq
