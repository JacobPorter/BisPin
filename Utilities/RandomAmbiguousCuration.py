#!/usr/bin/python
import optparse
import Constants
import datetime
import os
import sys

"""
@author: Jacob Porter

"""

def processACCFiles(directory_string):
    onlyfiles = [f for f in os.listdir(directory_string) if os.path.isfile(os.path.join(directory_string, f))]
    results_dictionary = {}
    counter = 0
    for f in onlyfiles:
        file_info = f.split(".")
        if file_info[-1] != "acc":
            continue
        i = int(file_info[-2])
        label = file_info[1]
        total_correct = extractStats(os.path.join(directory_string, f))
        if i in results_dictionary:
            results_dictionary[i][label] = total_correct
        else:
            counter += 1
            results_dictionary[i] = {label: total_correct}
    return (results_dictionary, counter)


def extractStats(filename_string):
    fd = open(filename_string)
    total = None
    correct = None
    for line in fd:
        if line.startswith("{"):
            for entry in line.split(","):
                if "Records_Analyzed" in entry:
                    total = entry.split(":")[1]
                elif "Correct" in entry:
                    correct = entry.split(":")[1]
    return (total, correct, int(correct) / (int(total) + 0.0))


def resultsString(results_dictionary):
    keys = results_dictionary.keys()
    keys.sort()
    results_string = ""
    first = True
    for k in keys:
        record_keys = results_dictionary[k].keys()            
        record_keys.sort()
        if first:
            first = False
            first_line = "Replicate\t"
            for rk in record_keys:
                first_line += rk + "_Total\t"
                first_line += rk + "_Correct\t"
                first_line += rk + "_Percent\t"
            results_string += first_line + "\n"
        line = str(k) + "\t"
        for rk in record_keys:
            line += "%s\t%s\t%f\t" % results_dictionary[k][rk]
        line += "\n"
        results_string += line
    return results_string


def main():
    now = datetime.datetime.now()
    usage = "usage: %prog [options] <Directory of ACC files from ambiguous simulations>"
    version = "%prog " + Constants.version
    description = ""
    epilog = Constants.creation_string
    p = optparse.OptionParser(usage = usage, version = version, 
                              description = description, epilog = epilog)
    _, args = p.parse_args()
    if len(args) != 1:
        p.error("There must be one argument.  Please check the input.")
    results_dictionary, counter = processACCFiles(args[0])
    sys.stdout.write(resultsString(results_dictionary))
    later = datetime.datetime.now()
    sys.stderr.write("Random Ambiguous Curation was run on the directory %s and found %d replicates.  The process took time %s.\n" % (args[0], counter, str(later - now)))

if __name__ == "__main__":
    main()