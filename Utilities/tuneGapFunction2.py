#!/usr/bin/python
import optparse
import os
import datetime
import subprocess
import sys
import multiprocessing
import shutil

"""
A file for processing failed BisPin with BFAST-Gap tuning because of missing SAM records.
"""

def deleteAllFiles(a_directory):
    if os.path.exists(a_directory):
        onlyfiles = [f for f in os.listdir(a_directory) if os.path.isfile(os.path.join(a_directory, f))]
        for f in onlyfiles:
            try:
                os.remove(os.path.join(a_directory, f))
            except OSError:
                pass
    return

def runBisPinCalculate(scoring_filename, num_threads, reference_genome_file, output_file, reads, total_reads, function_string, fstring_middle):
    outdir = os.path.dirname(output_file)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outdir_AUC = os.path.join(outdir, "AUC")
    outdir_source = os.path.join(outdir_AUC, "source")
    outdir_sink = os.path.join(outdir_AUC, "sink")
    if not os.path.exists(outdir_source):
        os.makedirs(outdir_source)
    if not os.path.exists(outdir_sink):
        os.makedirs(outdir_sink)
    num_processes = num_threads * 2
    report_files = [f for f in os.listdir(outdir) if os.path.isfile(os.path.join(outdir, f)) if ".BisPin.report" in f]
    if len(report_files) > 1:
        sys.stderr.write("A report file was found for:\t%s\n" % fstring_middle) #Won't show up for the subprocess.
        sys.stderr.flush()
        return fstring_middle
    else:
        sys.stderr.write("A report file was NOT found for:\t%s\n" % fstring_middle)
        sys.stderr.flush()
        #BisPin_arguments = ["BisPin_align.py", "-i", "2", "-x", scoring_filename, "-g", "-k", "-n", str(num_threads), "-b", "-1000", "-N", reference_genome_file, output_file, reads]
        #subprocess.call(BisPin_arguments, stdout = open(os.path.join(output_file + ".BisPin.out"), 'w'), stderr = open(os.path.join(output_file + ".BisPin.err"), 'w'))
        deleteAllFiles(outdir_AUC)
        deleteAllFiles(outdir_sink)
        onlyfiles = [f for f in os.listdir(outdir) if os.path.isfile(os.path.join(outdir, f)) if ".a4." in f and ".sam.gz" in f]
        for f in onlyfiles:
            shutil.copy2(os.path.join(outdir, f), os.path.join(outdir_source, f))
            subprocess.call(["gunzip", os.path.join(outdir_source, f)])
        #onlyfiles = [f for f in os.listdir(outdir) if os.path.isfile(os.path.join(outdir, f)) if ".sam" in f]
        #for f in onlyfiles:
        #    subprocess.Popen(["gzip", os.path.join(outdir, f)])
        calculateArguments = ["calculateAUCForBFAST_GAP.py", "-b", "0", "-e", "97", "-s", "1", "-o", outdir_sink, "-p", str(num_processes), reference_genome_file, outdir_source, reads, str(total_reads)]
        subprocess.call(calculateArguments, stdout = open(os.path.join(outdir_AUC, "AUC." + function_string + ".csv"), 'w'), stderr = open(os.path.join(outdir_AUC, "AUC." + function_string + ".err"), 'w' ))
        subprocess.Popen(["gzip", "-r", outdir_sink])
        if os.path.exists(outdir_source):
            deleteAllFiles(outdir_source)
            try:
                os.rmdir(outdir_source)
            except OSError:
                pass
        return fstring_middle
    
    
def parseArgs(options, args, p):
    if len(args) != 4:
        p.print_help()
        p.error("There must be four arguments.")
    extension = options.extension
    reference_genome_file = args[0]
    output_file = args[1]
    reads = args[2]
    number_of_reads = args[3]
    outdir = os.path.dirname(output_file)
    outfilename = os.path.basename(output_file)
    with open(options.scoring_file) as scoring_file:
        scoring_list = [n.strip() for n in scoring_file]
    fstring_begin = "_".join(scoring_list[0:5])
    if not extension:
        fstring_end = "_".join(scoring_list[9:len(scoring_list)])
    else:
        fstring_middle = "_".join(scoring_list[5:9])
    #gap_function_arguments = map(lambda x: float(x), scoring_list[5:9])
    gap_open_function_begin = [float(x) for x in options.gap_open_function_begin.split(",")]
    gap_open_function_end = [float(x) for x in options.gap_open_function_end.split(",")]
    gap_increment_function = [float(x) for x in options.gap_increment_function.split(",")]
    sys.stderr.write("Begin:\t%s\nEnd:\t%s\nStep:\t%s\n" % (str(gap_open_function_begin), str(gap_open_function_end), str(gap_increment_function)))
    sys.stderr.flush()
    pool = multiprocessing.Pool(processes = int(options.processes) / 2)
    processes_get = []
    i = gap_open_function_begin[0]
    counter = 0
    while i <= gap_open_function_end[0]:
        j = gap_open_function_begin[1]
        while j <= gap_open_function_end[1]:
            if not extension:
                fstring_middle = "%s_%s" % (str(i), str(j))
            else:
                fstring_end = "%s_%s" % (str(i), str(j))
            fstring = fstring_begin + "_" + fstring_middle + "_" + fstring_end
            my_output_file_name = outfilename + "." + fstring + ".sam"
            if not extension:
                fstring_write = fstring_middle
            else:
                fstring_write = fstring_end
            my_outfile = os.path.join(outdir, fstring_write, my_output_file_name)
            if not os.path.exists(os.path.join(outdir, fstring_write)):
                os.makedirs(os.path.join(outdir, fstring_write))
            scoring_filename = os.path.join(outdir, fstring_write, "scoring." + fstring + ".txt")
            with open(scoring_filename, 'w') as s_f:
                for k in range(5):
                    s_f.write(str(scoring_list[k]) + "\n")
                if not extension:
                    s_f.write(str(i) + "\n")
                    s_f.write(str(j) + "\n")
                    s_f.write(str(0.0) + "\n")
                    s_f.write(str(0.0) + "\n")
                    for k in range(4):
                        s_f.write(str(scoring_list[k + 9]) + "\n")
                    s_f.write(str(scoring_list[len(scoring_list) - 1]))
                else:
                    for k in range(5):
                        s_f.write(str(scoring_list[k + 5]) + "\n")
                    s_f.write(str(i) + "\n")
                    s_f.write(str(j) + "\n")
                    s_f.write(str(0.0) + "\n")
                    s_f.write(str(0.0) + "\n")
            p_handle = pool.apply_async(runBisPinCalculate, (scoring_filename, 2, reference_genome_file, my_outfile, reads, number_of_reads, fstring, fstring_write))
            sys.stderr.write("Added work to pool:\t%s\n" % (str(fstring_write)))
            sys.stderr.flush()
            processes_get.append(p_handle)
            counter += 1
            j += gap_increment_function[1]
        i += gap_increment_function[0]
    sys.stderr.write("Finished putting work on the pool.  There should be %d instances.\n" % (counter))
    sys.stderr.flush()
    completion_counter = 0
    for p in processes_get:
        completed_string = p.get()
        completion_counter += 1
        sys.stderr.write("Finished with work:\t%s\t%s\n" % (str(completed_string), str(completion_counter / (counter + 0.0))))
        sys.stderr.flush()
    
def main():
    now = datetime.datetime.now()
    usage = "usage: %prog [options] <reference_genome_file> <output_file> <reads> <number_of_reads>"
    description = ""
    epilog = ""
    p = optparse.OptionParser(usage = usage, 
                              description = description, epilog = epilog, )
    p.add_option('--scoring_file','-s',help='A scoring text file representing default values for the functions not being tuned. [default: %default]', default = "./scoring_function.txt")
    p.add_option('--gap_open_function_begin','-b',help='A string representing the beginning values of the function. [default: %default]', default = "0.05,-25.0")
    p.add_option('--gap_open_function_end','-e',help='A string representing the ending values of the function. [default: %default]', default = "1.25,150")
    p.add_option('--gap_increment_function', '-i', help='Amounts to increment each value. [default: %default]', default = "0.05,5")
    p.add_option('--processes', '-p', help = 'The number of processes to execute at once. [default: %default]', default = '8')
    p.add_option('--extension','-x',help='Tune the gap extension function instead of the gap open function.  [default: %default]', action='store_true', default = False)
    #p.add_option()
    options, args = p.parse_args()
    sys.stderr.write("Tuning on reads file:\t%s\n" % str(args[2]))
    sys.stderr.write("Starting at time:\t%s\n" % str(now))
    sys.stderr.flush()
    parseArgs(options, args, p)
    later = datetime.datetime.now()
    sys.stderr.write("Ending at time:\t%s\n" % str(later))
    sys.stderr.write("The tuning took time:\t%s\n" % (str(later - now)))

if __name__ == "__main__":
    main()
    
    