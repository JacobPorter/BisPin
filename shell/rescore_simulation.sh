#!/bin/sh
# $1 read iterate number
# $2 number of reads to simulate
~/Sherman_v0.1.7/Sherman -n $2  -e 2 -l 150 -CG 20 -CH 2 --genome_folder /research/jsporter/Data/genome/GRCm38.p5/genome_folder/ &> /research/jsporter/Data/reads/Mouse/Single/Simulation/simulated.rescore_test.single.$1.out

mv simulated.fastq /research/jsporter/Data/reads/Mouse/Single/Simulation/simulated.rescore_test.single.$1.fastq

BisPin_align.py -W -N -k -i 1 -n 6 -u /research/jsporter/Data/genome/GRCm38.p5/BisFAST/GRCm38.p5.fa /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.rescore_test.single.$1.sam /research/jsporter/Data/reads/Mouse/Single/Simulation/simulated.rescore_test.single.$1.fastq &> /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.rescore_test.single.$1.bispin.out

mkdir /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/rescore_test_$1

mv /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/*.sam /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/rescore_test_$1/

mv /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/rescore_test_$1/simulated.rescore_test.single.$1.sam /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.rescore_test.single.$1.sam

BisPin_postprocess.py -W /research/jsporter/Data/genome/GRCm38.p5/BisFAST/GRCm38.p5.fa /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/rescore_test_$1/ /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.rescore_test.rescored.single.$1.sam /research/jsporter/Data/reads/Mouse/Single/Simulation/simulated.rescore_test.single.$1.fastq

rm -r /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/rescore_test_$1

~/BisPin/Utilities/ExtractRescoredAlignments.py /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.rescore_test.rescored.single.$1.sam 1> /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.rescore_test.rescored_only.single.$1.sam

calculateSimulationAccuracy.py /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.rescore_test.rescored_only.single.$1.sam $2 1> /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.rescore_test.rescored_only.single.$1.acc 2> /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.rescore_test.rescored_only.single.$1.acc.out 

~/BisPin/Utilities/RandomAmbiguousChoice.py /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.rescore_test.single.$1.sam 1> /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.random_all_ambig.single.$1.sam

calculateSimulationAccuracy.py /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.random_all_ambig.single.$1.sam $2 1> /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.random_all_ambig.single.$1.acc 2> /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.random_all_ambig.single.$1.acc.out 

~/BisPin/Utilities/RandomAmbiguousChoice.py -f /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.rescore_test.rescored.single.$1.sam /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.rescore_test.single.$1.sam 1> /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.random_rescored_ambig.single.$1.sam

calculateSimulationAccuracy.py /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.random_rescored_ambig.single.$1.sam $2 1> /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.random_rescored_ambig.single.$1.acc 2> /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/simulated.random_rescored_ambig.single.$1.acc.out

mv /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/*.sam /research/jsporter/Data/output/Mouse/BisFAST/Single/Simulation/SAM/

#rm /research/jsporter/Data/reads/Mouse/Single/Simulation/simulated.single.$1.fastq
