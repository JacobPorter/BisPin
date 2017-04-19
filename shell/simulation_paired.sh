#!/bin/sh
# This produces the simulation file in the directory where the script is run
# $1 sample number
# $2 number of reads to simulate
# $3 Length of reads to simulate
~/Sherman_v0.1.7/Sherman -n $2 -e 2 -l $3 -CG 20 -CH 2 --paired_end --genome_folder /research/jsporter/Data/genome/GRCm38.p5/genome_folder/ &> /research/jsporter/Data/reads/Mouse/Paired-End/Simulation/simulated.paired.$1.$3.out

mv simulated_1.fastq /research/jsporter/Data/reads/Mouse/Paired-End/Simulation/simulated_1.paired.$1.$3.fastq

mv simulated_2.fastq /research/jsporter/Data/reads/Mouse/Paired-End/Simulation/simulated_2.paired.$1.$3.fastq

BisPin_align.py -i 1 -I 2,3 -n 6 -u /research/jsporter/Data/genome/GRCm38.p5/BisFAST/GRCm38.p5.fa /research/jsporter/Data/output/Mouse/BisFAST/Paired-End/Simulation/$3bp/simulated.paired.$1.sam /research/jsporter/Data/reads/Mouse/Paired-End/Simulation/simulated_1.paired.$1.$3.fastq /research/jsporter/Data/reads/Mouse/Paired-End/Simulation/simulated_2.paired.$1.$3.fastq   &> /research/jsporter/Data/output/Mouse/BisFAST/Paired-End/Simulation/$3bp/simulated.paired.$1.bispin.out

calculateSimulationAccuracy.py -p  /research/jsporter/Data/output/Mouse/BisFAST/Paired-End/Simulation/$3bp/simulated.paired.$1.sam $2 1> /research/jsporter/Data/output/Mouse/BisFAST/Paired-End/Simulation/$3bp/simulated.paired.$1.acc 2> /research/jsporter/Data/output/Mouse/BisFAST/Paired-End/Simulation/$3bp/simulated.paired.$1.acc.out

~/bismark_v0.16.3/bismark -p 6  --ambiguous --unmapped  --sam --temp_dir /research/jsporter/Data/output/Mouse/bismark/Paired-End/Simulation/ -o /research/jsporter/Data/output/Mouse/bismark/Paired-End/Simulation/$3bp/$1/    /research/jsporter/Data/genome/GRCm38.p5/bismark/ -1  /research/jsporter/Data/reads/Mouse/Paired-End/Simulation/simulated_1.paired.$1.$3.fastq -2  /research/jsporter/Data/reads/Mouse/Paired-End/Simulation/simulated_2.paired.$1.$3.fastq  &> /research/jsporter/Data/output/Mouse/bismark/Paired-End/Simulation/$3bp/simulated.paired.$1.bismark.out

calculateSimulationAccuracy.py -p  /research/jsporter/Data/output/Mouse/bismark/Paired-End/Simulation/$3bp/$1/simulated_1.paired.$1.$3_bismark_bt2_pe.sam $2 1> /research/jsporter/Data/output/Mouse/bismark/Paired-End/Simulation/$3bp/$1/simulated.paired.$1.acc 2> /research/jsporter/Data/output/Mouse/bismark/Paired-End/Simulation/$3bp/$1/simulated.paired.$1.acc.out

mv /research/jsporter/Data/output/Mouse/bismark/Paired-End/Simulation/$3bp/$1/simulated.paired.$1.acc /research/jsporter/Data/output/Mouse/bismark/Paired-End/Simulation/$3bp/simulated.paired.$1.acc

~/bwa-meth/bwameth.py --reference /research/jsporter/Data/genome/GRCm38.p5/bwa-meth/GRCm38.p5.fa /research/jsporter/Data/reads/Mouse/Paired-End/Simulation/simulated_1.paired.$1.$3.fastq /research/jsporter/Data/reads/Mouse/Paired-End/Simulation/simulated_2.paired.$1.$3.fastq  1> /research/jsporter/Data/output/Mouse/bwa-meth/Paired-End/Simulation/$3bp/simulated.paired.$1.sam 2> /research/jsporter/Data/output/Mouse/bwa-meth/Paired-End/Simulation/$3bp/simulated.paired.$1.bwameth.err

calculateSimulationAccuracy.py -p /research/jsporter/Data/output/Mouse/bwa-meth/Paired-End/Simulation/$3bp/simulated.paired.$1.sam $2 1> /research/jsporter/Data/output/Mouse/bwa-meth/Paired-End/Simulation/$3bp/simulated.paired.$1.acc 2> /research/jsporter/Data/output/Mouse/bwa-meth/Paired-End/Simulation/$3bp/simulated.paired.$1.acc.out

~/walt-1.0/bin/walt -t 6 -i /research/jsporter/Data/genome/GRCm38.p5/walt/GRCm38.p5.multiLine.fa.index.walt.dbindex -o /research/jsporter/Data/output/Mouse/walt/Paired-End/Simulation/$3bp/simulated.paired.$1.sam -1 /research/jsporter/Data/reads/Mouse/Paired-End/Simulation/simulated_1.paired.$1.$3.fastq  -2 /research/jsporter/Data/reads/Mouse/Paired-End/Simulation/simulated_2.paired.$1.$3.fastq  &> /research/jsporter/Data/output/Mouse/walt/Paired-End/Simulation/$3bp/simulated.paired.$1.walt.out

calculateSimulationAccuracy.py -p  /research/jsporter/Data/output/Mouse/walt/Paired-End/Simulation/$3bp/simulated.paired.$1.sam $2 1> /research/jsporter/Data/output/Mouse/walt/Paired-End/Simulation/$3bp/simulated.paired.$1.acc 2> /research/jsporter/Data/output/Mouse/walt/Paired-End/Simulation/$3bp/simulated.paired.$1.acc.out

rm /research/jsporter/Data/reads/Mouse/Paired-End/Simulation/simulated_1.paired.$1.$3.fastq

rm /research/jsporter/Data/reads/Mouse/Paired-End/Simulation/simulated_2.paired.$1.$3.fastq
