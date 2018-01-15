#!/bin/sh

~/BisPin/shell/trim_and_align_mouse.simulation.sh /research/jsporter/Data/reads/Mouse/Single/Simulation/simulated.rescore_test.single.500k.fastq /research/jsporter/Data/output/Mouse/trim/Simulation/simulated.sherman.500k 7 500000

~/BisPin/shell/trim_and_align_human.simulation.sh /research/jsporter/Data/reads/Human/IonTorrent/Simulation/reads.500k.dwgsim.human.ver0.fastq /research/jsporter/Data/output/Human/trim/Simulation/simulated.dwgsim.iontorrent.500k 7 500000
