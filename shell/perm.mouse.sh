#!/bin/sh

nohup ./perm10K /research/jsporter/Data/genome/GRCm38.p5/perm/GRCm38.p5.fa 101 --readFormat fastq --seed F3 -s /research/jsporter/Data/genome/GRCm38.p5/perm/GRCm38.p5.fa.index 1> /research/jsporter/Data/genome/GRCm38.p5/perm/perm.index.101bp.seedF3.out 2> /research/jsporter/Data/genome/GRCm38.p5/perm/perm.index.101bp.seedF3.err &

nohup ./perm10K /research/jsporter/Data/genome/GRCm38.p5/perm/GRCm38.p5.fa.index /research/jsporter/Data/reads/Mouse/Hairpin/BisFAST.hairpin.500k.R.fastq -o /research/jsporter/Data/output/Mouse/perm/recovered.perm.default.sam &> /research/jsporter/Data/output/Mouse/perm/recovered.perm.default.out &
