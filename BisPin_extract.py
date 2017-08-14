#!/usr/bin/python
import datetime
import sys
import optparse
from Utilities import IndentedHelpFormatterWithNL
from Utilities import SeqIterator
from Utilities import Constants

"""
This program extracts the uniquely mapped, ambiguously mapped, and unmapped SAM records into 
separate files.  The SAM file input is expected to be from a BFAST or BisPin SAM file. 
@author: Jacob Porter
"""

def isUnmapped(flag):
    """Checks if the SAM flag indicates an unmapped read."""
    return ((int(flag) >> 2) % 2) == 1

def isFiltered(flag):
    """Checks if the SAM flag indicates a filtered reads."""
    return ((int(flag) >> 9) % 2) == 1

def extract(sam_file, output_file, no_output):
    sam_input = SeqIterator.SeqIterator(sam_file, file_type=Constants.SAM)
    if not no_output:
        output_unique_name = output_file + ".unique.sam"
        output_ambig_name = output_file + ".ambig.sam"
        output_unmap_name = output_file + ".unmap.sam"
        output_filt_name = output_file + ".filt.sam"
        output_unique = SeqIterator.SeqWriter(open(output_unique_name, 'w'), file_type=Constants.SAM)
        output_ambig = SeqIterator.SeqWriter(open(output_ambig_name, 'w'), file_type=Constants.SAM)
        output_unmap = SeqIterator.SeqWriter(open(output_unmap_name, 'w'), file_type=Constants.SAM)
        output_filt = SeqIterator.SeqWriter(open(output_filt_name, 'w'), file_type=Constants.SAM)
    else:
        output_unique = None
        output_ambig = None
        output_unmap = None
        output_filt = None
    last_record = []
    counts = {"unique" : 0, "unmap": 0, "ambig": 0, "filt": 0}
    for record in sam_input:
        if isUnmapped(record[Constants.SAM_KEY_FLAG]):
            last_record = writeTo(last_record, output_unique, output_ambig, counts, no_output)
            counts["unmap"] += 1
            if not no_output:
                output_unmap.write(record)
        elif isFiltered(record[Constants.SAM_KEY_FLAG]):
            last_record = writeTo(last_record, output_unique, output_ambig, counts, no_output)
            counts["filt"] += 1
            if not no_output:
                output_filt.write(record)
        elif len(last_record) >= 1 and record[Constants.SAM_KEY_QNAME] != last_record[0][Constants.SAM_KEY_QNAME]:
            last_record = writeTo(last_record, output_unique, output_ambig, counts, no_output)
            last_record.append(record)
        else:
            last_record.append(record)
    writeTo(last_record, output_unique, output_ambig, counts, no_output)
    return counts
                    
def writeTo(last_record, output_unique, output_ambig, counts, no_output):
    if len(last_record) == 1:
        if not no_output:
            output_unique.write(last_record[0])
        counts["unique"] += 1
    elif len(last_record) > 1:
        if not no_output:
            for ambig_record in last_record:
                        output_ambig.write(ambig_record)
        counts["ambig"] += 1
    return []
            
            
def main():
    """
    Divides a BFAST or BisPin SAM file into three separate SAM files: unique, ambiguous, unmapped.
    The percentage of each type is returned.
    """
    now = datetime.datetime.now()
    usage = "usage: %prog [options] <sam_file> <output_file> "
    version = "%prog " + Constants.version
    description = ""
    epilog = ""
    p = optparse.OptionParser(usage = usage, version = version, 
                              description = description, epilog = epilog, 
                              formatter = IndentedHelpFormatterWithNL.IndentedHelpFormatterWithNL())
    p.add_option('--no_output', '-N', help='Do not create the output files. [default: %default]', action='store_true', default=False)
    options, args = p.parse_args()
    if len(args) < 1:
        p.print_help()
        p.error("There must be at least one argument.")
    elif options.no_output:
        output_file = None
    elif len(args) != 2:
        p.error("There must be two arguments.")
    else:
        output_file = args[1]
    sam_file = args[0]
    sys.stderr.write("Extracting reads for file %s and outputting results to files starting with %s\n" % (sam_file, str(output_file)))
    counts = extract(sam_file, output_file, options.no_output)
    total = 0
    for key in counts:
        total += counts[key]
    if total == 0:
        total = 0.000000000000000000000000000000001
    else:
        total += 0.0
    sys.stderr.write("Unique:\t%d\t%f\nAmbig:\t%d\t%f\nUnmap:\t%d\t%f\nFilt:\t%d\t%f\nTotal:\t%d\n" % (counts["unique"], 100 * counts["unique"] / total, counts["ambig"], 100 * counts["ambig"] / total, counts["unmap"], 100 * counts["unmap"] / total, counts["filt"], 100 * counts["filt"] / total, total))
    later = datetime.datetime.now()
    runtime = later - now
    sys.stderr.write("Runtime:\t" + str(runtime) + "\n")

if __name__ == "__main__":
    main()