#!/usr/bin/python
import extractInformation
import os
import sys
import optparse
import datetime
import csv

"""
@author: Jacob Porter

nohup createInfoTrimSheet.py . 1> /research/jsporter/Data/output/Human/InfoTrim/result_compilation/InfoTrim.simulation.results.tsv 2> /research/jsporter/Data/output/Human/InfoTrim/result_compilation/InfoTrim.simulation.results.err &
"""

def getArguments(file_name, programs = ["bispin", "bismark", "bwameth", "walt"]):
    file_name_lower = file_name.lower()
    program_name = None
    for p in programs:
        if p in file_name_lower:
            program_name = p
            break
    file_name_list = file_name.replace("_", ".").strip().split(".")
    r = -1.0
    s = -1.0
    m = -1
    #sys.stderr.write("file_name_list:\t" + str(file_name_list) + "\n")
    for i in range(len(file_name_list)):
        if file_name_list[i].startswith("r-"):
            r = float(file_name_list[i].split("-")[1] + "." + file_name_list[i + 1])
        elif file_name_list[i].startswith("s-"):
            s = float(file_name_list[i].split("-")[1] + "." + file_name_list[i + 1])
        elif file_name_list[i].startswith("m-"):
            m = int(file_name_list[i].split("-")[1])
    return (program_name, r, s, m)

def processResultFiles(directory_string, argFunc):
    onlyfiles = [f for f in os.listdir(directory_string) if os.path.isfile(os.path.join(directory_string, f))]
    accfiles = [f for f in onlyfiles if f.endswith(".acc")]
    resultsfiles = [f for f in onlyfiles if not f.endswith(".acc") and (f.endswith(".txt") or f.endswith(".report") or f.endswith(".count") or f.endswith(".mapstats"))]
    accfiles_dictionary = {argFunc(acc) : extractInformation.extractAccuracy(os.path.join(directory_string, acc)) for acc in accfiles}
    results_dictionary = { }
    for f in resultsfiles:
        if f.endswith(".txt") and "bismark" in f: #bismark
            results_dictionary[argFunc(f)] = extractInformation.extractBismark(os.path.join(directory_string, f))
        elif f.endswith(".report"): #BisPin
            results_dictionary[argFunc(f)] = extractInformation.extractBisPin(os.path.join(directory_string, f))
        elif f.endswith(".count"): #BWAMeth
            results_dictionary[argFunc(f)] = extractInformation.extractBWAMeth(os.path.join(directory_string, f))
        elif f.endswith(".mapstats"): #Walt
            results_dictionary[argFunc(f)] = extractInformation.extractWalt(os.path.join(directory_string, f))
    return (accfiles_dictionary, results_dictionary)

def createSheet(accfiles_dictionary, results_dictionary):
    fieldnames = ["Program", "r", "s", "m", "Reads", "Unique", "Ambig", "Unmap", "Filt", "Unmap + Filt", "Total", "Recall", "Precision", "F1_Score", "CpG", "CHG", "CHH"]
    sheetwriter = csv.DictWriter(sys.stdout, fieldnames, restval='', extrasaction='ignore', delimiter='\t')
    sheetwriter.writeheader()
    for key in results_dictionary:
        write_dictionary = results_dictionary[key]
        write_dictionary["Program"] = key[0]
        write_dictionary["r"] = key[1]
        write_dictionary["s"] = key[2]
        write_dictionary["m"] = key[3]
        if key in accfiles_dictionary:
            write_dictionary.update(accfiles_dictionary[key])
        else:
            sys.stderr.write("The key %s was not found in the accuracy files dictionary, but it was present in the results dictionary.\n" % (str(key)))
        sheetwriter.writerow(write_dictionary)
        
def main():
    now = datetime.datetime.now()
    usage = "usage: %prog [options] <Directory of InfoTrim simulation results>"
    description = ""
    p = optparse.OptionParser(usage = usage, 
                              description = description)
    _, args = p.parse_args()
    if len(args) != 1:
        p.error("There must be one argument.  Please check the input.")
    if not os.path.isdir(args[0]):
        p.error("The directory could not be found.")
    sys.stderr.write("Extracting information from files located in the directory %s." % str(args[0]))
    accfiles_dictionary, results_dictionary = processResultFiles(args[0], getArguments)
    createSheet(accfiles_dictionary, results_dictionary)
    later = datetime.datetime.now()
    sys.stderr.write("The process took time %s.\n" % (str(later - now)))

if __name__ == "__main__":
    main()