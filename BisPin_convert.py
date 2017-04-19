#!/usr/bin/python
from Utilities import BisPin_util
from Utilities import Constants
from Utilities import IndentedHelpFormatterWithNL
import optparse
import os
import sys
import multiprocessing
import datetime

"""
@author: Jacob Porter
@summary: A program for converting bisulfite DNA read data where C's are converted to T's or G's are converted to A's
"""

logstr = "\nBisPin_convert: "

converter_description="\
BisPin_convert creates two copies of the \
supplied DNA sequence file.  They are stored in  \
the same directory as the supplied file.  \
One copy is C-to-T converted and the other copy is G-to-A converted.\n\n\
Required:\n\
\t is the DNA sequence file to convert.  It can be either a FASTQ or FASTA file.\n\
"

def main():
    """
    The entry point into the program.
    """
    usage = "usage: %prog [options] <DNA_sequence_file>"
    version = "%prog " + Constants.version
    p = optparse.OptionParser(usage = usage, 
                              version = version, 
                              description = converter_description,
                              epilog = Constants.creation_string,
                              formatter = IndentedHelpFormatterWithNL.IndentedHelpFormatterWithNL())
    p.add_option('--type' , '-t', help='The type of file to process. [default: %default]\nOptions:\nFASTA\nFASTQ', default=Constants.FASTA)
    options, args = p.parse_args()
    if len(args) == 0:
        p.print_help()
    if len(args) != 1:
        p.error("There is one required argument.")
    fastafile = args[0]
    filename = os.path.basename(fastafile)
    directory = os.path.dirname(fastafile)
    file_type = options.type.lower()
    if file_type != Constants.FASTA and file_type != Constants.FASTQ:
        p.error("The file type is unrecognized.  Please check.")
    if not os.path.exists(fastafile):
        p.error("%sThe file at %s could not be located." % (logstr, fastafile))
    sys.stderr.write("%sProcessing %s sequence file at %s\n" % (logstr, file_type, fastafile))
    sys.stderr.flush()
    now = datetime.datetime.now()
    p_CT = multiprocessing.Process(target=BisPin_util.convert_seqs, args=(directory, filename, True, None, file_type))
    p_GA = multiprocessing.Process(target=BisPin_util.convert_seqs, args=(directory, filename, False, None, file_type))
    p_CT.start()
    p_GA.start()
    p_CT.join()
    finish_string = "%sFinished creating %s sequence file at location %s.\n"
    sys.stderr.write(finish_string % (logstr, "C to T", directory))
    p_GA.join()
    sys.stderr.write(finish_string % (logstr, "G to A", directory))
    sys.stderr.flush()
    later = datetime.datetime.now()
    sys.stderr.write("%sElapsed time -- %s" % (logstr, str(later - now)))
    sys.stderr.flush()
    
    
if __name__ == '__main__':
    main()