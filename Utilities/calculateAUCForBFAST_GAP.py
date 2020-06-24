#!/usr/bin/python
import datetime
import optparse
import multiprocessing
import subprocess
import sys
import os
import csv
from . import calculateSimulationAccuracy

"""
This program examines all filter values in some range and computes the F1_Score based for simulated reads.  It uses the preprocessed BFAST SAM files.
@author: Jacob Porter
"""

def callPostProcess(filter_value, ref_genome, sam_input_directory, reads_file, output_directory, useSecondaryIndexes, total_reads, deleteMe, sherman):
    """
    Creates a postprocessed SAM file and computes F1_Score, precision, and recall for a given filter value cutoff.
    """
    filename = os.path.basename(reads_file)
    output_file = os.path.join(output_directory, filename + "." + filter_value + ".sam")
    postprocess_file_string_out = os.path.join(output_directory, filename + "." + filter_value + ".out")
    postprocess_file_string_err = os.path.join(output_directory, filename + "." + filter_value + ".err")
#     sys.stderr.write("Filename:\t%s\nOutput_File:\t%s\nSAM_Input:\t%s\n" % (filename, output_file, sam_input_directory))
#     sys.stderr.flush()
    postprocess_args = ["BisPin_postprocess.py",  "-b",  filter_value, "-N"]
    if useSecondaryIndexes:
        postprocess_args.append("-i")
    postprocess_args += [ref_genome, sam_input_directory, output_file, reads_file]
    subprocess.call(postprocess_args, stdout = open(postprocess_file_string_out, 'w'), stderr = open(postprocess_file_string_err, 'w'))
    counter_dict = calculateSimulationAccuracy.processSAM(output_file, False, 3, False, not sherman, print_value = 0)
    reads_analyzed = int(counter_dict["Records_Analyzed"]) + 0.0
    if reads_analyzed == 0.0:
        reads_analyzed = 0.0000000000000001
    recall = counter_dict["Correct"] / total_reads
    prec = counter_dict["Correct"] / reads_analyzed
    if prec <= 0.00000000000001 or recall <= 0.00000000000001:
        F1_score = 0.0
    else:
        F1_score = 2* prec * recall / (prec + recall)
    if deleteMe:
        os.remove(output_file)
        os.remove(postprocess_file_string_out)
        os.remove(postprocess_file_string_err)
    return (filter_value, F1_score, prec, recall)

def parseArgs(options, args):
    """
    Parses the options and args and dispatches processes to postprocess the BFAST SAM files.
    Collects the F1_Score, AUC, etc. and writes the results to a CSV file to stdout.
    """
    pool = multiprocessing.Pool(processes = int(options.processes))
    begin = int(options.begin)
    end = int(options.end)
    step = int(options.step)
    ref_genome = args[0]
    sam_input_directory = args[1]
    reads_file = args[2]
    processes_get = []
    sherman = options.sherman
    for i in range(begin, end, step):
        p_handle = pool.apply_async(callPostProcess, (str(i), ref_genome, sam_input_directory, reads_file, options.output_directory, options.useSecondaryIndexes, float(args[3]), not options.no_delete, sherman))
        processes_get.append(p_handle)
    max_F1_score = -1.0
    max_filter_value = None
    AUC = 0.0
    fieldnames = ["Filter_Value", "F1_Score", "Precision", "Recall", "AUC", "Max_F1_Score", "Maximizing_Filter_Value"]
    sheetwriter = csv.writer(sys.stdout, delimiter='\t')
    sheetwriter.writerow(fieldnames)
    counter = 0
    for p in processes_get:
        counter += 1
        info = p.get()
        (filter_value, F1_score, _, _) = info
        if F1_score > max_F1_score:
            max_F1_score = F1_score
            max_filter_value = filter_value
        AUC += F1_score * step
        sheetwriter.writerow(info)
        sys.stdout.flush()
    end_row = ["", "", "", ""]
    sheetwriter.writerow(end_row + [AUC, max_F1_score, max_filter_value])
    return (AUC, max_F1_score, max_filter_value, counter)
        
def getReport(results, timing, initial_string, useInitial):
    """
    Returns a string giving a summary and report of the execution, AUC, etc.
    """
    AUC, max_F1_score, max_filter_value, counter = results
    report = "AUC:\t%s\nmax_F1_score:\t%s\nmax_filter_value:\t%s\ncounter:\t%s\ntiming:\t%s\nJacob Porter Virginia Tech\n" % (str(AUC), str(max_F1_score), str(max_filter_value), str(counter), str(timing))
    if useInitial:
        report = initial_string + report
    return report
    
def main():
    """
    Gives the command line interface, checks the arguments and options, dispatches the parser, and prints the report to stderr.
    """
    now = datetime.datetime.now()
    usage = "usage: %prog [options] <reference_genome_file> <sam_input_directory> <reads_file> <total_reads> "
    description = "This program examines all filter values in some range and computes the F1_Score based for simulated reads.  It uses the preprocessed BFAST SAM files."
    epilog = ""
    p = optparse.OptionParser(usage = usage, 
                              description = description, epilog = epilog)
    p.add_option('--processes', '-p', help='The number of processes to use to get results. [default: %default]', default = '12')
    p.add_option('--begin', '-b', help='The beginning filter value. [default: %default]', default = '0')
    p.add_option('--end', '-e', help='The ending filter value. [default: %default]', default = '97')
    p.add_option('--step', '-s', help='The step increment for the filter value. [default: %default]', default = '1')
    p.add_option('--output_directory', '-o', help='The directory to store temporary files and write results. [default: %default]', default = './')
    p.add_option('--useSecondaryIndexes','-i',help='Turn this flag on if the SAM files were created with secondary indexes.  [default: %default]', action='store_true', default = False)
    p.add_option('--no_delete','-n',help='Do not delete the postprocessed SAM files and BisPin output.  Deletes the files by default. [default: %default]', action='store_true', default = False)
    p.add_option('--sherman','-m',help='The reads files were generated by sherman (rather than DWGSIM). [default: %default]', action='store_true', default = False)
    options, args = p.parse_args()
    if len(args) != 4:
        p.print_help()
        p.error("There must be four arguments.")
    if not os.path.isdir(options.output_directory) or not os.path.isdir(args[1]):
        p.error("The output directory or the input direcotry does not exist.")
    if not os.path.exists(args[0]) or not os.path.exists(args[2]):
        p.error("The reference genome file or the reads file does not exist.")
    try:
        float(args[3])
        int(options.processes)
        float(options.begin)
        float(options.end)
        float(options.step)
    except ValueError:
        p.error("A value that should be a numeric constant is not.")
    initial_string = "The process calculateAUCForBFAST_GAP was started on:\t%s\nReference genome:\t%s\nThe input:\t%s\nThe reads file:\t%s\nTotal reads:\t%s\n" % (str(now), args[0], args[1], args[2], args[3])
    sys.stderr.write(initial_string)
    sys.stderr.flush()
    results = parseArgs(options, args)
    later = datetime.datetime.now()
    report = getReport(results, later - now, initial_string, False)
    sys.stderr.write("%s" % (report))
    
if __name__ == "__main__":
    main()