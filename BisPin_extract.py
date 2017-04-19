#!/usr/bin/python
import datetime
import sys
import optparse
from Utilities import IndentedHelpFormatterWithNL
from Utilities import SeqIterator
from Utilities import Constants

"""
@author: Jacob Porter
TODO: Does this work?
TODO: Improve the sequence iterator to extract BisPin fields.
"""

def extract(sam_file, output_file):
    sam_input = SeqIterator.SeqIterator(sam_file, file_type=Constants.SAM)
    output_unique_name = output_file + ".unique.sam"
    output_ambig_name = output_file + ".ambig.sam"
    output_unmap_name = output_file + ".unmap.sam"
    output_unique = SeqIterator.SeqWriter(open(output_unique_name, 'w'), file_type=Constants.SAM)
    output_ambig = SeqIterator.SeqWriter(open(output_ambig_name, 'w'), file_type=Constants.SAM)
    output_unmap = SeqIterator.SeqWriter(open(output_unmap_name, 'w'), file_type=Constants.SAM)
    last_record = []
    for record in sam_input:
        if record[Constants.SAM_KEY_FLAG] == Constants.SAM_VALUE_UNMAPPED:
            last_record = writeTo(last_record, output_unique, output_ambig)
            output_unmap.write(record)
        elif len(last_record) >= 1 and record[Constants.SAM_KEY_QNAME] != last_record[0][Constants.SAM_KEY_QNAME]:
            last_record = writeTo(last_record, output_unique, output_ambig)
            last_record.append(record)
        else:
            last_record.append(record)
    writeTo(last_record, output_unique, output_ambig)
                    
def writeTo(last_record, output_unique, output_ambig):
    if len(last_record) == 1:
        output_unique.write(last_record[0])
    elif len(last_record) > 1:
        for ambig_record in last_record:
                    output_ambig.write(ambig_record)
    return []
            
            
def main():
    """
    Divides a BFAST or BisPin SAM file into three separate SAM files: unique, ambiguous, unmapped.
    """
    now = datetime.datetime.now()
    usage = "usage: %prog [options] <sam_file> <output_file> "
    version = "%prog " + Constants.version
    description = ""
    epilog = ""
    p = optparse.OptionParser(usage = usage, version = version, 
                              description = description, epilog = epilog, 
                              formatter = IndentedHelpFormatterWithNL.IndentedHelpFormatterWithNL())
    _, args = p.parse_args()
    if len(args) != 2:
        p.print_help()
        p.error("There must be two arguments.")
    sam_file = args[0]
    output_file = args[1]
    sys.stderr.write("Extracting reads for file %s and outputting results to files starting with %s\n" % (args[0], args[1]))
    extract(sam_file, output_file)
    later = datetime.datetime.now()
    runtime = later - now
    sys.stderr.write(str(runtime) + "\n")

if __name__ == "__main__":
    main()