#!/usr/bin/python
from . import SeqIterator
import optparse
import datetime
import os
import sys

"""
@author: Jacob Porter

"""

def selectRecords(readsfile, begin, end, file_type):
    reads_iterator = SeqIterator.SeqIterator(readsfile, file_type = file_type)
    output_writer = SeqIterator.SeqWriter(sys.stdout, file_type = file_type)
    counter = 0
    for record in reads_iterator:
        if file_type == 'fastq' or file_type == 'fasta':
            record_length = len(record[1])
        else:
            record_length = len(record["SEQ"])
        if record_length >= begin and record_length <= end:
            output_writer.write(record)
            counter += 1
    output_writer.flush()
    output_writer.close()
    return counter


def postInfo(readsfile, file_type, count, now, later):
    sys.stderr.write("reads file:\t%s\nfile type:\t%s\ncount:\t%s\nStart and end times:\t%s\t%s\n" % (readsfile, file_type, str(count), str(now), str(later)))


def main():
    now = datetime.datetime.now()
    usage = "usage: %prog [options] <min_length> <max_length>  <input_file_path> "
    version = "%prog "
    description = ""
    epilog = ""
    p = optparse.OptionParser(usage = usage, version = version, 
                              description = description, epilog = epilog)
    p.add_option('--file_type','-t',help='The file type.  Either FASTA or FASTQ. [default: %default]', default = "SAM")
    options, args = p.parse_args()
    if len(args) != 3:
        p.error("There are not enough arguments.  Please check the input.")
    readsfile = args[2]
    if not os.path.exists(readsfile):
        p.error("The reads file does not exist or could not be accessed.")
    file_type = options.file_type.lower()
    if file_type != "fastq" and file_type != "fasta" and file_type != "sam":
        p.error("The file type must be either FASTA, FASTQ, or SAM.")
    try:
        begin = int(args[0])
        end = int(args[1])
    except TypeError:
        p.error("The minimum and maximum lengths must be integers.")
    except ValueError:
        p.error("The minimum and maximum lengths must be integers.")
    if not begin <= end or not begin >= 0 or not end >= 0:
        p.error("Both the mimimum and maximum lengths must be larger than zero, and the minimum length must be less than or equal to the ending length.")
    count = selectRecords(readsfile, begin, end, file_type)
    later = datetime.datetime.now()
    postInfo(readsfile, file_type, count, now, later)


if __name__ == "__main__":
    main()