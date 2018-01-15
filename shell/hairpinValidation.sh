#!/bin/sh
# $1 test sam file

PREFIX="hairpin.R.1mil"

hairpinValidation.py $PREFIX.bfast.sam $1 1> $1.validate.bfast.val 2> $1.validate.bfast.out &
hairpinValidation.py $PREFIX.bowtie2.sam $1 1> $1.validate.bowtie2.val 2> $1.validate.bowtie2.out &
hairpinValidation.py $PREFIX.bwa.sam $1 1> $1.validate.bwa.val 2> $1.validate.bwa.out &
