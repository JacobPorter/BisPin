#!/bin/sh
# $1
/home/jsporter/bfast-gap-0.1.0/bfast/bfast match -f /research/jsporter/Data/genome/GRCh38.p9/BisFAST/GRCh38.p9.fa.BisPin.CT -i 1 -r /research/jsporter/Data/reads/Human/RegularIonTorrent/test.fastq -T /research/jsporter/Data/output/Human/test/ 1> /research/jsporter/Data/output/Human/test/test.bfast-gap.bmf 2> /research/jsporter/Data/output/Human/test/test.bfast-gap.err

/home/jsporter/bfast-gap-0.1.0/bfast/bfast localalign -f /research/jsporter/Data/genome/GRCh38.p9/BisFAST/GRCh38.p9.fa.BisPin.CT -m /research/jsporter/Data/output/Human/test/test.bfast-gap.bmf -x ~/bfast-gap-0.1.0/scoring_function.txt  1> /research/jsporter/Data/output/Human/test/test.bfast-gap.baf 2> /research/jsporter/Data/output/Human/test/test.bfast-gap.align.err

/home/jsporter/bfast-gap-0.1.0/bfast/bfast postprocess -f /research/jsporter/Data/genome/GRCh38.p9/BisFAST/GRCh38.p9.fa.BisPin.CT -i /research/jsporter/Data/output/Human/BFAST-Gap/test/test.bfast-gap.baf -x ~/bfast-gap-0.1.0/scoring_function.txt 1> /research/jsporter/Data/output/Human/BFAST-Gap/test/test.bfast-gap.sam 2> /research/jsporter/Data/output/Human/BFAST-Gap/test/test.bfast-gap.post.err
