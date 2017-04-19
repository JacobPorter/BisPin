#!/usr/bin/python
import datetime
import sys
import optparse
import random
import IndentedHelpFormatterWithNL
import SeqIterator
import Constants

"""
@author: Jacob Porter
"""

def extract(sam_file):
    sam_input = SeqIterator.SeqIterator(sam_file, file_type=Constants.SAM)
    sam_writer = SeqIterator.SeqWriter(sys.stdout, file_type=Constants.SAM)
    last_record = []
    rescored_count = 0
    for record in sam_input:
        if record[Constants.SAM_KEY_FLAG] == Constants.SAM_VALUE_UNMAPPED:
            continue
        elif len(last_record) >= 1 and record[Constants.SAM_KEY_QNAME] != last_record[0][Constants.SAM_KEY_QNAME]:
            rescored_count += writeTo(last_record, sam_writer)
            last_record = [record]
        else:
            last_record.append(record)
    rescored_count += writeTo(last_record, sam_writer)
    sam_writer.flush()
    return rescored_count


def writeTo(last_record, sam_writer):
    if len(last_record) == 1 and Constants.SAM_KEY_RESCORE in last_record[0]:
        sam_writer.write(last_record[0])
        return 1 
    else:
        return 0
            
            
def main():
    """
    """
    now = datetime.datetime.now()
    usage = "usage: %prog [options] <sam_BisPin_file> "
    version = "%prog " + Constants.version
    description = ""
    epilog = ""
    p = optparse.OptionParser(usage = usage, version = version, 
                              description = description, epilog = epilog, 
                              formatter = IndentedHelpFormatterWithNL.IndentedHelpFormatterWithNL())
    options, args = p.parse_args()
    if len(args) != 1:
        p.print_help()
        p.error("There must be a SAM file of ambiguous reads.")
    sam_file = args[0]
    sys.stderr.write("")
    rescored_count = extract(sam_file)
    later = datetime.datetime.now()
    runtime = later - now
    sys.stderr.write("Rescored count:\t%s\nRuntime:\t%s\n" % (str(rescored_count), str(runtime)))

if __name__ == "__main__":
    main()