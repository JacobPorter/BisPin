# BisPin
BisPin is a Python 2.7 program that uses BFAST to map bisulfite-treated short DNA reads to a reference genome.  It supports the hairpin construction strategy and rescores ambiguously mapped reads.  Its default alignment scoring function is biologically motivated.

Created by Jacob Porter while at Virginia Tech.
http://www.jacobporter.com
http://www.vt.edu

## Requirements
BisPin requires BFAST, which can be gotten at https://sourceforge.net/projects/bfast/.
BisPin requires Python 2.7, which can be gotten at https://www.python.org/downloads/.

## Using BisPin
If BFAST is not in the PATH variable, then it needs to be specified as an option to the BisPin index and align programs.  See the online help documentation for each program for usage and options.  This can be accessed by typing, for example, "BisPin_align.py --help".  The workflow for BisPin is the following.

1. BisPin_convert.py  This program creates two copies of the reference genome.  One is C-to-T converted, and the other is G-to-A converted.  This is necessary for the bisulfite alignment.  The files are stored in the same directory as the reference genome file.
2. BisPin_index.py  This program creates BFAST indexes for the converted genomes in the previous step.  If the hairpin recovery strategy is used, BFAST indexes for the regular genome need to be created with this program as well.  The indexes are stored in the same directory as the reference genome.  The converted genomes from the previous step must exist and be in the same directory as the reference genome.
3. BisPin_aign.py This program calls BFAST to do the alignments on the reads file.  The output is a single SAM file that contains the uniquely mapped, ambiguously mapped, and unmapped reads, and a summary report is written as well.  For more information on the SAM file format see https://samtools.github.io/hts-specs/.

## Extra processing options
If desired, the BisPin\_align.py program will save the intermediate BFAST SAM files with the -k option.  The BisPin\_postprocess.py program can be used to experiment with different postprocessing options by taking these BFAST SAM files as input without incurring the cost of alignment each time.

The program BisPin_extract.py can be used to separate the uniquely mapped, ambiguously mapped, and unmapped partitions of the SAM file.

##Known Issues
- The support for multiple index functionality when using a partition and multiprocessing is unimplemented.
- Support for input file partitioning for the hairpin recovery processing is not implemented.
- The hairpin recovered methylation calling is inconsistent with compared to not using hairpin recovery.
- Add comments and test BisPin_extract
- BFAST gives incorrect flags for paired end mapping that indicate incorrect orientations.
            #These flags are for the GA converted genome
            #115 should be 83 RC
            #179 should be 147 not RC
            #These flags are for the CT converted genome
            #65 should be 99 not RC
            #129 should be 147 RC