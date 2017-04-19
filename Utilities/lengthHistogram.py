#!/usr/bin/python
import SeqIterator
import optparse
import datetime
import os
import sys
import math

"""
@author: Jacob Porter
This program calculates a histogram vector for read lengths.
"""

def calculateHistogram(readsfile, file_type):
    reads_iterator = SeqIterator.SeqIterator(readsfile, file_type = file_type)
    histogram = {}
    counter = 0
    for record in reads_iterator:
        counter += 1
        record_length = len(record[1])
        if record_length in histogram:
            histogram[record_length] += 1
        else:
            histogram[record_length] = 1
    keys = histogram.keys()
    keys.sort()
    items = []
    average = 0
    median = -1
    number_so_far = 0
    sum_of_squares = 0
    for k in keys:
        items.append(histogram[k])
        average += histogram[k] * k
        sum_of_squares += k * k * histogram[k]
        number_so_far += histogram[k]
        if median == -1 and number_so_far >= (counter / 2.0):
            median = k
    std_dev = math.sqrt((counter * sum_of_squares - average * average) / (counter * (counter - 1)))
    average = average / (counter + 0.0)
    sys.stdout.write(str(keys) + "\n")
    sys.stdout.write(str(items) + "\n")
    sys.stdout.write("Average:\t%s\n" % (str(average)))
    sys.stdout.write("Median:\t%s\n" % (str(median)))
    sys.stdout.write("Standard Deviation:\t%s\n" % (str(std_dev)))
    return counter


def postInfo(readsfile, file_type, count, now, later):
    sys.stderr.write("%s\t%s\t%s\t%s\t%s\t\n" % (readsfile, file_type, str(count), str(now), str(later)))


def main():
    now = datetime.datetime.now()
    usage = "usage: %prog [options] <reads_file_path> "
    version = "%prog "
    description = ""
    epilog = ""
    p = optparse.OptionParser(usage = usage, version = version, 
                              description = description, epilog = epilog)
    p.add_option('--file_type','-t',help='The file type.  Either FASTA or FASTQ. [default: %default]', default = "FASTQ")
    options, args = p.parse_args()
    if len(args) != 1:
        p.error("A file path must be specified.  Please check the input.")
    readsfile = args[0]
    if not os.path.exists(readsfile):
        p.error("The reads file does not exist or could not be accessed.")
    file_type = options.file_type.lower()
    if file_type != "fastq" and file_type != "fasta":
        p.error("The file type must be either FASTA or FASTQ.")
    count = calculateHistogram(readsfile, file_type)
    later = datetime.datetime.now()
    postInfo(readsfile, file_type, count, now, later)


if __name__ == "__main__":
    main()