#!/usr/bin/python
import datetime
import sys
import optparse
import random
from . import IndentedHelpFormatterWithNL
from . import SeqIterator
from . import Constants

"""
@author: Jacob Porter
"""

def runSimulation(norescore_sam_file, rescore_sam_file):
    norescore_ambig_dict = extract(norescore_sam_file)
    rescore_ambig_dict = None
    if rescore_sam_file != None:
        rescore_ambig_dict = extract(rescore_sam_file)
    simulate(norescore_ambig_dict, rescore_ambig_dict)


def extract(sam_file):
    sam_input = SeqIterator.SeqIterator(sam_file, file_type=Constants.SAM)
    last_record = []
    ambig_dictionary = {}
    for record in sam_input:
        if record[Constants.SAM_KEY_FLAG] == Constants.SAM_VALUE_UNMAPPED:
            continue
        elif len(last_record) >= 1 and record[Constants.SAM_KEY_QNAME] != last_record[0][Constants.SAM_KEY_QNAME]:
            writeTo(last_record, ambig_dictionary)
            last_record = [record]
        else:
            last_record.append(record)
    writeTo(last_record, ambig_dictionary)
    return ambig_dictionary


def writeTo(last_record, ambig_dictionary):
    if len(last_record) > 1:
        ambig_dictionary[last_record[0][Constants.SAM_KEY_QNAME]] = last_record
    

def simulate(norescore_ambig_dict, rescore_ambig_dict):
    SAM_writer = SeqIterator.SeqWriter(sys.stdout, file_type=Constants.SAM)
    for QNAME in norescore_ambig_dict:
        if rescore_ambig_dict == None or QNAME not in rescore_ambig_dict:
            record_list = norescore_ambig_dict[QNAME]
            selected_record = record_list[random.randint(1, len(record_list)) - 1]
            SAM_writer.write(selected_record)
    SAM_writer.flush()
           
            
def main():
    """
    """
    now = datetime.datetime.now()
    usage = "usage: %prog [options] <sam_ambiguous_file> "
    version = "%prog " + Constants.version
    description = ""
    epilog = ""
    p = optparse.OptionParser(usage = usage, version = version, 
                              description = description, epilog = epilog, 
                              formatter = IndentedHelpFormatterWithNL.IndentedHelpFormatterWithNL())
    p.add_option('--AmbigFile', '-f', help='The SAM file produced by BisPin where rescoring has occured.', default=None)
    options, args = p.parse_args()
    if len(args) != 1:
        p.print_help()
        p.error("There must be a SAM file of ambiguous reads.")
    sam_file = args[0]
    sys.stderr.write("")
    runSimulation(sam_file, options.AmbigFile)
    later = datetime.datetime.now()
    runtime = later - now
    sys.stderr.write(str(runtime) + "\n")

if __name__ == "__main__":
    main()