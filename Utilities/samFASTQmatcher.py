#!/usr/bin/python
import optparse
import datetime
import sys
import SeqIterator

"""
This program takes a SAM file and a FASTQ file as input.  It prints out fastq file records
where the fastq id matches some QNAME in the SAM file.
@author: Jacob Porter
"""

def samMatchFASTQ(sam_file, fasta_file):
    sam_dictionary = {record["QNAME"] : True for record in SeqIterator.SeqIterator(sam_file, file_type='SAM')}
    fasta_iterator = SeqIterator.SeqIterator(fasta_file, file_type = 'fastq')
    fasta_writer = SeqIterator.SeqWriter(sys.stdout, file_type = 'fastq')
    counter = 0
    #sys.stderr.write("sam_dictionary:\n%s\n" % str(sam_dictionary))
    for fasta_record in fasta_iterator:
#         sys.stderr.write("FASTQ record:\t%s\n" % fasta_record[0])
#         sys.stderr.flush()
        if sam_dictionary.get(fasta_record[0], False):
            counter += 1
            fasta_writer.write(fasta_record)
    return counter

def parseArgs(p, now):
    _, args = p.parse_args()
    if len(args) != 2:
        p.print_help()
        p.error("There must be at least two arguments.")
    sys.stderr.write("Finding FASTQ records from %s that match the QNAME from SAM records from %s\n" % (args[0], args[1]))
    sys.stderr.flush()
    counter = samMatchFASTQ(args[0], args[1])
    sys.stderr.write("Found %d FASTQ records that matched a record in the SAM file.\n" % counter)
    later = datetime.datetime.now()
    sys.stderr.write("The process took time:\t%s\n" % str(later - now))
    sys.stderr.flush()
    
    
def main():
    now = datetime.datetime.now()
    usage = "usage: %prog [options] <sam_file> <fastq_file> "
    description = "This program takes a SAM file and a FASTQ file as input.  It prints out fastq file records where the fastq id matches some QNAME in the SAM file."
    epilog = ""
    p = optparse.OptionParser(usage = usage, 
                              description = description, epilog = epilog)
    parseArgs(p, now)
    #p.add_option('--processes', '-p', help='The number of processes to use to get results. [default: %default]', default = '12')
    
if __name__ == "__main__":
    main()