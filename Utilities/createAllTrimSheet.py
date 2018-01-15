#!/usr/bin/python
import createInfoTrimSheet
import extractInformation
import os
import sys
import optparse
import datetime
import csv
import subprocess
import copy
import multiprocessing

"""
For all the trimmers, extract mapping quality information.
@author: Jacob Porter
"""

def getTrimArguments(file_name, trimmers = ["notrim", "trimmomatic", "sickle", "reaper", "infotrim", "erne", "cutadapt"], mappers = ["bispin", "bismark", "bwameth", "walt"]):
    my_trimmer = None
    my_mapper = None
    for trim in trimmers:
        if trim in file_name:
            my_trimmer = trim
            break
    for prog in mappers:
        if prog in file_name:
            my_mapper = prog
            break
    if my_trimmer == None:
        my_trimmer = "notrim"
    return (my_trimmer, my_mapper)

def createSheet(accfiles_dictionary, results_dictionary, hairpin_dictionary):
    fieldnames = ["Trimmer", "Mapper", "Reads", "Unique", "Ambig", "Unmap", "Filt", "Unmap + Filt", "Total", "Recall", "Precision", "F1_Score", "CpG", "CHG", "CHH", "Hairpin_Hits_bwa", "Test_Correct_bwa", "Test_Incorrect_bwa", "Hairpin_Hits_bfast", "Test_Correct_bfast", "Test_Incorrect_bfast", "Hairpin_Hits_bowtie2", "Test_Correct_bowtie2", "Test_Incorrect_bowtie2"]
    sheetwriter = csv.DictWriter(sys.stdout, fieldnames, restval='', extrasaction='ignore', delimiter='\t')
    sheetwriter.writeheader()
    for key in results_dictionary:
        write_dictionary = results_dictionary[key]
        write_dictionary["Trimmer"] = key[0]
        write_dictionary["Mapper"] = key[1]
        if key in accfiles_dictionary:
            write_dictionary.update(accfiles_dictionary[key])
        else:
            sys.stderr.write("The key %s was not found in the accuracy files dictionary, but it was present in the results dictionary.\n" % (str(key)))
        if key in hairpin_dictionary:
            write_dictionary.update(hairpin_dictionary[key])
        else:
            sys.stderr.write("The key %s was not found in the hairpin validation files dictionary, but it was present in the results dictionary.\n" % (str(key)))
        sheetwriter.writerow(write_dictionary)

def callSimAcc(directory_string, myArgs, sfile):
        subprocess.call(myArgs, stdout = open(os.path.join(directory_string, sfile + ".acc"), 'w'), stderr = open(os.path.join(directory_string, sfile + ".arr"), 'w'))

def getValFilesDictionary(directory_string, argFunc):
    #def defaultValDictionary():
    #    return {"Hairpin_Hits_bwa": 0, "Hairpin_Hits_bfast": 0, "Hairpin_Hits_bowtie2": 0, "Test_Correct_bwa" : 0, "Test_Correct_bfast" : 0, "Test_Correct_bowtie2" : 0, "Test_Incorrect_bwa" : 0, "Test_Incorrect_bfast" : 0, "Test_Incorrect_bowtie2" : 0}
    def defaultDict():
        return {"Test_Correct" : 0.0, "Test_Incorrect": 0.0, "counter": 0.0, "F1_score": 0.0}
    regular_mappers = [".bwa.", ".bfast.", ".bowtie2."]
    val_fields_to_extract = ["Hairpin_Hits", "Test_Correct", "Test_Incorrect"]
    onlyfiles = [f for f in os.listdir(directory_string) if os.path.isfile(os.path.join(directory_string, f))]
    valfiles = [f for f in onlyfiles if f.endswith(".val")]
    if valfiles == []:
        return {}
    val_dictionary = {argFunc(f) : {} for f in valfiles}
    avg_dictionary = {"notrim" : defaultDict(), "trimmomatic": defaultDict(), "sickle" : defaultDict(), "reaper" : defaultDict(), "infotrim": defaultDict(), "erne" : defaultDict(), "cutadapt": defaultDict()}
    for f in valfiles:
        my_rm = None
        for rm in regular_mappers:
            if rm in f:
                my_rm = rm[1:len(rm) - 1]
                break
        extracted_dictionary = extractInformation.extractHairpin(f)
        key = argFunc(f)
        my_trimmer = key[0]
        entry_dictionary = val_dictionary[argFunc(f)]
        for val_field in val_fields_to_extract:
            selector = 1
            if val_field.startswith("Hairpin_Hits"):
                selector = 0
            else:
                avg_dictionary[my_trimmer][val_field] += float(extracted_dictionary[val_field][selector])
                avg_dictionary[my_trimmer]["counter"] += 1.0
            entry_dictionary[val_field + "_" + my_rm] = extracted_dictionary[val_field][selector]
    for trimmer in avg_dictionary:
        avg_dictionary[trimmer]["counter"] = (avg_dictionary[trimmer]["counter"] / 2.0)
        avg_dictionary[trimmer]["Test_Correct"] /= avg_dictionary[trimmer]["counter"]
        avg_dictionary[trimmer]["Test_Incorrect"] /= avg_dictionary[trimmer]["counter"]
    return val_dictionary, avg_dictionary

def getAllACCFiles(directory_string, num_reads, dwgsim, pool):
    onlyfiles = [f for f in os.listdir(directory_string) if os.path.isfile(os.path.join(directory_string, f))]
    samfiles = [f for f in onlyfiles if f.endswith(".sam") and ".a4." not in f]
    simAccArgs = ["calculateSimulationAccuracy.py"]
    if dwgsim:
        simAccArgs.append("-d")
    pd_list = []
    for sfile in samfiles:
        myArgs = copy.copy(simAccArgs)
        myArgs +=  [os.path.join(directory_string, sfile), str(num_reads)]
        pd = pool.apply_async(callSimAcc, (directory_string, myArgs, sfile))
        pd_list.append(pd)
    for pd in pd_list:
        pd.get()
        
def write_avg(avg_dictionary):
    fieldnames = ["Trimmer", "Test_Correct" , "Test_Incorrect", "counter", "F1_score"]
    sheetwriter = csv.DictWriter(open("hairpin.averages.csv", 'w'), fieldnames, restval='', extrasaction='ignore', delimiter='\t')
    sheetwriter.writeheader()
    for key in avg_dictionary:
        write_dictionary = avg_dictionary[key]
        write_dictionary["Trimmer"] = key
        sheetwriter.writerow(write_dictionary)

def main():
    now = datetime.datetime.now()
    usage = "usage: %prog [options] <Directory of trimming results>"
    description = ""
    p = optparse.OptionParser(usage = usage, 
                              description = description)
    
    p.add_option('--num_reads', '-N', help='The total number of reads from the input file.  This is for computing accuracy results. [default: %default]', default = 500000)
    p.add_option('--computeAccuracy', '-a', help='Compute the accuracy files.  Use this when the accuracy files do not exist.', action='store_true', default=False)
    p.add_option('--dwgsim', '-d', help='This flag should be set to true if the simulations were generated with dwgsim and the accuracy files need to be computed.', action='store_true', default=False)
    options, args = p.parse_args()
    if len(args) != 1:
        p.error("There must be one argument.  Please check the input.")
    if not os.path.isdir(args[0]):
        p.error("The directory could not be found.")
    sys.stderr.write("Extracting information from files located in the directory %s." % str(args[0]))
    if options.computeAccuracy:
        getAllACCFiles(args[0], options.num_reads, options.dwgsim, multiprocessing.Pool(processes = 12))
    accfiles_dictionary, results_dictionary = createInfoTrimSheet.processResultFiles(args[0], getTrimArguments)
    hairpin_dictionary, avg_dictionary = getValFilesDictionary(args[0], getTrimArguments)
    write_avg(avg_dictionary)
    createSheet(accfiles_dictionary, results_dictionary, hairpin_dictionary)
    later = datetime.datetime.now()
    sys.stderr.write("The process took time %s.\n" % (str(later - now)))

if __name__ == "__main__":
    main()