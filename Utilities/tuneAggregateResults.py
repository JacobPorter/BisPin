#!/usr/bin/python
import optparse
import os
import datetime
import sys
import csv

"""
This gets the information from multiple tuning executions and writes a csv file to stdout summarizing the information.
@author: Jacob Porter
"""

def extractInfoFromFile(file_location, fstring):
    output_dictionary = {"F_string": fstring, "AUC": None, "max_F1_score": None, "max_filter_value": None, "Total reads": None}
    with open(file_location) as my_file:
        for line in my_file:
            for key in output_dictionary:
                if line.startswith(key):
                    output_dictionary[key] = float(line.split(":")[1])
    return output_dictionary
    
def extractInfoFromDirectory(root_directory):
    sub_directories = next(os.walk(root_directory))[1]
    info_dictionary = {}
    counter = 0
    max_AUC = -1
    maximizers_AUC = []
    max_F1 = -1
    maximizers_F1 = []
    for sd in sub_directories:
        file_directory = os.path.join(root_directory, sd, "AUC")
        if os.path.isdir(file_directory):
            info_file = [f for f in os.listdir(file_directory) if os.path.isfile(os.path.join(file_directory, f)) and f[::-1].startswith("rre.")][0]
            info_dictionary[sd] = extractInfoFromFile(os.path.join(file_directory, info_file), sd)
            if float(info_dictionary[sd]["AUC"]) > max_AUC:
                max_AUC = float(info_dictionary[sd]["AUC"])
                maximizers_AUC = [sd]
            elif float(info_dictionary[sd]["AUC"]) == max_AUC:
                maximizers_AUC.append(sd)
            if float(info_dictionary[sd]["max_F1_score"]) > max_F1:
                max_F1 = float(info_dictionary[sd]["max_F1_score"])
                maximizers_F1 = [sd]
            elif float(info_dictionary[sd]["max_F1_score"]) == max_F1:
                maximizers_F1.append(sd)
            counter += 1
    sys.stderr.write("Processed %d files.\n" % counter)
    sys.stderr.write("The max AUC was:\t%s\n" % str(max_AUC))
    sys.stderr.write("The maximizing AUC functions were:\t%s\n" % ",".join(maximizers_AUC))
    sys.stderr.write("The max F1 score was:\t%s\n" % str(max_F1))
    sys.stderr.write("The functions that maximized the F1 score were:\t%s\n" % ",".join(maximizers_F1))
    sys.stderr.flush()
    return info_dictionary

def writeTuningInfo(info_dictionary):
    with sys.stdout as out:
        tuning_writer = csv.DictWriter(out, fieldnames = ["F_string", "AUC", "max_F1_score", "max_filter_value", "Total reads"], delimiter = "\t", quotechar = '"')
        tuning_writer.writeheader()
        for fstring in info_dictionary:
            tuning_writer.writerow(info_dictionary[fstring])

def parseArgs(options, args, p):
    if len(args) != 1:
        p.print_help()
        p.error("There must be one argument, the directory.")
    if not os.path.isdir(args[0]):
        p.error("The root directory does not exists.  Please check the argument.")
    sys.stderr.write("Starting to aggregate tuning information for directory:\t%s\n" % args[0])
    sys.stderr.flush()
    writeTuningInfo(extractInfoFromDirectory(args[0]))

def main():
    now = datetime.datetime.now()
    usage = "usage: %prog [options] <root_directory>"
    description = "This gets the information from multiple tuning executions and writes a csv file to stdout."
    epilog = ""
    p = optparse.OptionParser(usage = usage, 
                              description = description, epilog = epilog, )
    options, args = p.parse_args()
    parseArgs(options, args, p)
    later = datetime.datetime.now()
    sys.stderr.write("The process took time:\t%s\n" % (str(later - now)))
    sys.stderr.flush()

if __name__ == "__main__":
    main()
    
    
