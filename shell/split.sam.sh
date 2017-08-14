#!/bin/sh
# $1 SAM file
# $2 prefix for outputfile
# $3 begin length
# $4 end length
lengthSelector.py $3 $4 $1 1> $2.$3-$4.sam 2> $2.$3-$4.lengthSelector.err

BisPin_extract.py -N $2.$3-$4.sam 1> $2.$3-$4.count 2> $2.$3-$4.count.err

alignmentScoreDistribution.py -b 20 $2.$3-$4.sam 1> $2.$3-$4.scores 2> $2.$3-$4.scores.err
