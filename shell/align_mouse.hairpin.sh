#!/bin/sh
# $1 The reads directory and filename 
# $2 The output directory and filename
# $3 The number of threads

~/bwa-meth/bwameth.py --reference /research/jsporter/Data/genome/GRCm38.p5/bwa-meth/GRCm38.p5.fa $1 1> $2.bwameth.sam 2> $2.bwameth.err

countBWA.py $2.bwameth.sam 1> $2.bwameth.count 2> $2.bwameth.count.out

~/BisPin/shell/hairpinValidation.sh $2.bwameth.sam

~/walt-1.0/bin/walt -t $3 -i /research/jsporter/Data/genome/GRCm38.p5/walt/GRCm38.p5.multiLine.fa.index.walt.dbindex -o $2.walt.sam -r $1  &> $2.walt.out

~/BisPin/shell/hairpinValidation.sh $2.walt.sam

BisPin_align.py -n $3 -k -i 2 /research/jsporter/Data/genome/GRCm38.p5/BisFAST/GRCm38.p5.fa $2.bispin.sam $1 &> $2.bispin.out

~/BisPin/shell/hairpinValidation.sh $2.bispin.sam

~/bismark_v0.16.3/bismark -p $3  --sam -o $2.bismark.sam  /research/jsporter/Data/genome/GRCm38.p5/bismark/ $1  &> $2.bismark.out

~/BisPin/shell/hairpinValidation.sh $2.bismark.sam

