#!/bin/sh

# $1 input fastq file
# $2 reference genome file
# $3 scoring matrix file
# $4 number of threads
# $5 output prefix

bfast match -T ./ -n $4 -f $2 -i 2 -r $1 1> $5.bmf 2> $5.bmf.out
bfast localalign -n $4 -f $2 -x $3 -m $5.bmf 1> $5.baf 2> $5.baf.out
bfast postprocess -O 1 -n $4 -f $2 -x $3 -a 4  -i $5.baf  1> $5.sam 2> $5.sam.out
