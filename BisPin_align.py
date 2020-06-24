#!/usr/bin/python
import sys
sys.dont_write_bytecode = True
import optparse
import os
import multiprocessing
import subprocess
import math
import itertools
import string
import random
import datetime
import BisPin_postprocess
from Utilities import SeqIterator
from Utilities import SeqDoubleIterator
from Utilities import BisPin_util
from Utilities import IndentedHelpFormatterWithNL
from Utilities import Constants
#from Utilities import Deprecator

"""
@author: Jacob Porter
@summary: Aligns bisulfite-treated reads from a FASTQ file to a reference genome FASTA file with BFAST.
@requires:  The indices must already be created. 
TODO: add support for input file partitioning to the hairpin recovery processing with multiple indexes and parallel postprocessing.
TODO: For BFAST-Gap, the initialization values for the Smith-Waterman matrix are not using the gap open / extension run length weights but are using the default values.
"""

logstr = "\nBisPin_align: "

def run_align(path_to_bfast, reads_file, fasta_file, fasta_dir, tmpDir, numThreads, type_str, scoring_file, layout, additional_parameters, usePipe, recover=False):
    """
    Calls BFAST to align the FASTQ DNA sequences to the reference genome.
    This function is meant to run in parallel in its own process.
    A sempaphore is used to allow running this function sequentially rather than in parallel.
    @param path_to_bfast: A string giving the executable for BFAST
    @param reads_file: A string giving the location of the FASTQ reads file
    @param fasta_file: The file name only for the reference genome.
    @param fasta_dir: The directory where the reference genome is located.
    @param tmpDir: The location of a temporary directory where itermediate and temporary BFAST files will be written too.
    @param numThreads: The number of threads to run the alignment process for BFAST
    @param type_str: A string giving the conversion type of the reads and the reference genome.  See Constants for examples.
    @param scoring_file: The location of the file that determines the alignment scoring function.
    @param additional_parameters: A dictionary that stores additional parameters to control BFAST.
    @param recover: Used for hairpin data.  This controls whether the recovery step will be used.
    @return: A list of temporary files to delete and timing information for the BFAST execution phases.
    """
    #Semaphore for controlling how many processes are running at a time.
    global alignment_semaphore
    alignment_semaphore.acquire(block = True, timeout = None)
    #Extract extra arguments
    gzip_str = "-z" if additional_parameters.get(Constants.GZIP, False) else ""
    insertSizeAvg = additional_parameters.get(Constants.insertSizeAvg, None)
    insertSizeStdDev = additional_parameters.get(Constants.insertSizeStdDev, None)
    keySize = additional_parameters[Constants.keySize]
    maxKeyMatches = additional_parameters[Constants.maxKeyMatches]
    maxNumReadMatches = additional_parameters[Constants.maxNumReadMatches]
    maxNumAlignmentMatches = additional_parameters[Constants.maxNumAlignmentMatches]
    mainIndexes = additional_parameters[Constants.mainIndexes]
    secondaryIndexes = additional_parameters[Constants.secondaryIndexes]
    offsets = additional_parameters[Constants.offsets]
    hairpin_paired = additional_parameters[Constants.HAIRPIN_PAIRED]
    #Get the correct reference genome FASTA file
    (read_str, fasta_str) = BisPin_util.getReadStrFastaStr(type_str)
    if recover:
        fasta_path = os.path.join(fasta_dir, fasta_file)
    else:
        fasta_path = os.path.join(fasta_dir, fasta_file + ".BisPin." + fasta_str)
    read_path = reads_file
    reads_file = os.path.basename(reads_file)
    #Construct the BFAST matches arguments.
    matches_filename = os.path.join(tmpDir, reads_file + "." + type_str + '.matches.bmf')
    if (type_str.startswith(Constants.CONV_CT) and type_str.endswith(Constants.CONV_CT)) or (type_str.startswith(Constants.CONV_GA) and type_str.endswith(Constants.CONV_GA)):
        direction_of_alignments = 1
    else:
        direction_of_alignments = 2
    if recover:
        direction_of_alignments = 1  #Is this correct?
    matches_params = [path_to_bfast, 'match', gzip_str, '-f', fasta_path, '-r', read_path, '-w', str(direction_of_alignments), '-T', tmpDir, '-n', numThreads, '-t']
    if keySize != None:
        matches_params.append('-k')
        matches_params.append(str(keySize))
    if maxKeyMatches != None:
        matches_params.append('-K')
        matches_params.append(str(maxKeyMatches))
    if maxNumReadMatches != None:
        matches_params.append('-M')
        matches_params.append(str(maxNumReadMatches))
    if mainIndexes != None:
        matches_params.append('-i')
        matches_params.append(str(mainIndexes))
    if secondaryIndexes != None:
        matches_params.append('-I')
        matches_params.append(str(secondaryIndexes))
    if offsets != None:
        matches_params.append('-o')
        matches_params.append(str(offsets))
    #Construct the BFAST localalign arguments.
    align_filename = os.path.join(tmpDir, reads_file + "." + type_str  + '.aligned')
    localalign_args = [path_to_bfast, 'localalign', '-f', fasta_path, '-x', scoring_file, '-n', numThreads, '-t' , '-m', matches_filename]
    if maxNumAlignmentMatches != None:
        localalign_args.append('-M')
        localalign_args.append(str(maxNumAlignmentMatches))
    #Find matches and do alignments with BFAST.  Write the output to a file.
    sys.stderr.write("%sCalling BFAST on %s threads for the %s reads and %s reference.\n" % (logstr, numThreads, read_str, fasta_str))
    if not usePipe:
        sys.stderr.write("%sFinding matches for %s reads.\n" % (logstr, read_str))
        sys.stderr.flush()
        match_now = datetime.datetime.now()
        matches_fd = open(matches_filename, 'w')
        subprocess.call(matches_params, stdout = matches_fd)
        match_later = datetime.datetime.now()
        match_elapsed = match_later - match_now
        matches_fd.close()
        #Do the alignments with BFAST localalign
        sys.stderr.write("%sFinding alignments for %s reads.\n" % (logstr, read_str))
        sys.stderr.flush()
        align_now = datetime.datetime.now()
        align_fd = open(align_filename, 'w')
        subprocess.call(localalign_args, stdout = align_fd)
        align_fd.close()
        align_later = datetime.datetime.now()
        align_elapsed = align_later - align_now
    else: #Pipe matches to localalign
        align_now = datetime.datetime.now()
        localalign_args = localalign_args[0 : len(localalign_args) - 2]
        align_fd = open(align_filename, 'w')
        sys.stderr.write("%s BFAST match and alignment parameters: %s\t%s\n" % (logstr, matches_params, localalign_args))
        match_popen = subprocess.Popen(matches_params, stdout = subprocess.PIPE)
        align_popen = subprocess.Popen(localalign_args, stdin = match_popen.stdout, stdout = align_fd)
        align_popen.communicate()
        align_fd.close()
        align_later = datetime.datetime.now()
        align_elapsed = align_later - align_now
    #Generate a raw BFAST SAM file with postprocess
    sys.stderr.write("%sPerforming BFAST postprocessing for %s reads.\n" % (logstr, read_str))
    sys.stderr.flush()
    post_now = datetime.datetime.now()
    if layout == Constants.LAYOUT_SINGLE or recover or hairpin_paired == False:
        postprocess_args = [path_to_bfast, 'postprocess', '-f', fasta_path, '-i', align_filename, '-a', '4', '-x', scoring_file, '-Y', '2', '-n', numThreads, '-t'] 
    else:
        postprocess_args = [path_to_bfast, 'postprocess', '-f', fasta_path, '-i', align_filename, '-a', '4', '-x', scoring_file, '-Y', '0', '-P', '0', '-S', '0', '-n', numThreads, '-t']
    if insertSizeAvg != None:
        postprocess_args.append("-v")
        postprocess_args.append(str(float(insertSizeAvg)))
    if insertSizeStdDev != None:
        postprocess_args.append("-s")
        postprocess_args.append(str(float(insertSizeStdDev)))
    postprocess_file_string = samPath(tmpDir, reads_file + ".a4", read_str, fasta_str)
    subprocess.call(postprocess_args, stdout = open(postprocess_file_string, 'w'))
    post_later = datetime.datetime.now()
    post_elapsed = post_later - post_now
    if not usePipe:
        timing_info = {"bfast_post_processing": post_elapsed, "bfast_match": match_elapsed, "bfast_alignment": align_elapsed}
    else:
        timing_info = {"bfast_post_processing": post_elapsed, "bfast_match_and_alignment": align_elapsed}
    sys.stderr.write("%sFinished BFAST postprocessing for %s reads.\n" % (logstr, read_str))
    sys.stderr.flush()
    #Release the semaphore and return temporary file names and timing information
    alignment_semaphore.release()
    if not usePipe:
        return ([postprocess_file_string, matches_filename, align_filename], timing_info)
    else:
        return ([postprocess_file_string, align_filename], timing_info)


def init_child(semaphore_):
    """
    This is for the process pool for alignments.
    It initializes the semaphore used to run groups of processes at a time.
    @param: semaphore_: a multiprocessing semaphore object.
    """
    global alignment_semaphore
    alignment_semaphore = semaphore_


def makeOneFilePairedEnd(reads1, reads2, output_file_name, gzip_switch, noRC):
    """
    Since BFAST requires a single file, this function converts two files representing paired sequence data into a single file.
    @param reads1: The forward reads FASTQ file name.
    @param reads2: The reverse reads FASTQ file name.  These will be reverse complemented.
    @param output_file_name: A string that indicates the location of the FASTQ file to write the results to.
    @param gzip_switch: If True, then the file represented by output_file_name will be gzip compressed.  If False, no compression will be used.
    @param noRC: When set to True, the reverse complement of the reads2 file will NOT be taken.  If False, then the reverse complement will be taken.  
    """
    paired_ends = SeqDoubleIterator.SeqDoubleIterator(reads1, reads2, file_type=Constants.FASTQ, gzip_switch = gzip_switch)
    output_writer = SeqIterator.SeqWriter(BisPin_util.getPossibleGZIPFile(output_file_name, gzip_switch), 
                                          file_type=Constants.FASTQ)
    for record in paired_ends:
        seq_id = BisPin_util.extractPairedID(record[0][0], record[1][0])
        record1 = list(record[0])
        record1[0] = seq_id
        record2 = list(record[1])
        record2[0] = seq_id
        if not noRC:
            record2[1] = BisPin_util.reverseComplement(record2[1])
            record2[2] = record2[2][::-1]
        output_writer.write(record1)
        output_writer.write(record2)


def samPath(tmpDir, reads_file, read_str, fasta_str):
    """
    Gets the path to a SAM file.  This is used to write the SAM file output by BFAST.
    @param tmpDir: the temporary directory where the SAM file is
    @param reads_file: the file name for a read file.  This will be the basis of the SAM file name.
    @param read_str: a string representing how the read was converted.
    @param fasta_str: a string representing how the reference genome was converted.
    @return: A string representing a path to a SAM file.
    """
    out_str = read_str + "_" + fasta_str + ".sam"
    return os.path.join(tmpDir, reads_file + "." + out_str)


def createAlignmentScoringFunctionFile(file_directory, gap_open, gap_extension, nucleotide_match, nucleotide_mismatch):
    """
    Creates a file that represents the affine scoring function for scoring alignments.
    @param file_directory: A string representing the directory to store the scoring file.
    @param gap_open:  The affine gap open penalty.
    @param gap_extension: The affine gap extension penalty
    @param nucleotide_match: The nucleotide match score.
    @param nucleotide_mismatch:  The nucleotide mismatch penalty.    
    @return: A string representing the file path.
    """
    rand_string = "".join([random.choice(string.ascii_letters + string.digits) for _ in range(Constants.BISPIN_RANDOM_LENGTH)])
    file_name = Constants.ALIGNFILENAME + rand_string + '.txt'
    file_path = os.path.join(file_directory, file_name)
    file_object = open(file_path, 'w')
    file_object.write(str(gap_open) + "\n")
    file_object.write(str(gap_extension) + "\n")
    file_object.write(str(nucleotide_match) + "\n")
    file_object.write(str(nucleotide_mismatch) + "\n")
    file_object.close()
    return file_path
    

def recover(CtoT_file, GtoA_file, tmpDir, gzip_switch, start, end):
    """
    Performs the hairpin recovery strategy.
    @param CtoT_file: The C to T converted hairpin file.
    @param GtoA_file: The G to A converted hairpin file.
    @param tmpDir: The directory to write the resulting files to.
    @return: Returns a list of temporary files created.
    """
    hairpin_rand_string = "".join([random.choice(string.ascii_letters) for _ in range(Constants.BISPIN_RANDOM_LENGTH)])
    outFile = os.path.join(tmpDir, Constants.HAIRPINFILEPREFIX + hairpin_rand_string)
    CtoT_Seqs = SeqIterator.SeqIterator(CtoT_file, file_type = Constants.FASTQ, gzip_switch = gzip_switch)
    GtoA_Seqs = SeqIterator.SeqIterator(GtoA_file, file_type = Constants.FASTQ, gzip_switch = gzip_switch)
    sys.stderr.write("%sRecovering original sequences from the T-enriched file %s and the A-enriched file %s and writing output to %s\n" % (logstr, CtoT_file, GtoA_file, outFile))
    sys.stderr.flush()
    R_filename = os.path.join(outFile + '.R') 
    T_filename = os.path.join(outFile + '.T')
    A_filename = os.path.join(outFile + '.A')
    post_filename = os.path.join(outFile + '.post')
    fd_R = BisPin_util.getPossibleGZIPFile(R_filename, gzip_switch)
    fd_T = BisPin_util.getPossibleGZIPFile(T_filename, gzip_switch)
    fd_A = BisPin_util.getPossibleGZIPFile(A_filename, gzip_switch)
    orig_writer = SeqIterator.SeqWriter(fd_R, file_type = Constants.FASTQ)
    T_writer = SeqIterator.SeqWriter(fd_T, file_type = Constants.FASTQ)
    A_writer = SeqIterator.SeqWriter(fd_A, file_type = Constants.FASTQ)
    fd_post = open(post_filename, 'w')
    counter = -1
    for i,j in zip(CtoT_Seqs, GtoA_Seqs):
        counter += 1
        if (start != None and counter < start):
            continue
        elif (end != None and counter > end):
            break
        if sameSeq(i,j):
            align(i, j, orig_writer, T_writer, A_writer, fd_post)
        else:
            sys.stderr.write("%sError!!! Sequences do not have the same identifier.  %s, %s\n" % (logstr, i[0], j[0]))
            sys.stderr.flush()
    return (R_filename, T_filename, A_filename, post_filename)


def align(i, j, orig_writer, T_writer, A_writer, fd_post):
    """
    This function determines if paired hairpin sequences match up perfectly or not.
    If they do, then the recovered sequence is written out, and if not, the two sequences
    are written out.
    @param i: One sequence record (a tubple for the FASTQ file)
    @param j: Another sequence record.
    @param fd_R: file object for writing the recovered sequences.
    @param fd_T: file object for writing the C to T converted sequences that do not recover
    @param fd_A: file object for writing the G to A converted sequences that do not recover
    @return: (my_align, mismatches, mismatch_cnt, percent_match)
    The tuple gives the alignment, a string representing the location and identity of mismatches,
    the number of mismatches, and the percent of the sequence that matches.
    """
    len1 = len(i[1])
    len2 = len(j[1])
    my_align = ""
    mismatches = ""
    mismatch_cnt = 0
    for pos in range(min(len1, len2)):
        ib = i[1][pos].upper()
        jb = j[1][pos].upper()
        if ib == jb:
            my_align += ib
        elif ib == "T" and jb == "C":
            my_align += "C"
        elif ib == "G" and jb == "A":
            my_align += "G"
        elif ib != jb:
            my_align += "*"
            mismatches += str(pos) + "_" + ib + jb + ":"
            mismatch_cnt += 1
    percent_match = (max(len1, len2) - mismatch_cnt - math.fabs(len1-len2)) / (max(len1, len2) + 0.0)
    if percent_match == 1.0:
        orig_writer.write((i[0], my_align, BisPin_util.averagePhredStrings(i[2][0:len(my_align)], j[2][0:len(my_align)])))
        fd_post.write("%s\t%s\n" % (str(i[0]), Constants.CONV_R_R))
    else:
        T_writer.write(i)
        A_writer.write(j)
        fd_post.write("%s\t%s\n" % (str(i[0]), Constants.CONV_CT))
    return (my_align, mismatches, mismatch_cnt, percent_match)
        

def sameSeq(i,j):
    '''
    For the hairpin data.  This determines if the ids for two sequences are identical.
    @param i: One sequence record
    @param j: Another sequence record
    @return: A boolean indicating whether the sequences are identical.
    '''
    return i[0] == j[0]


def deleteTempFiles(list_of_temp_files, deleteMe = True):
    """
    Deletes the temporary files created by the alignment process.
    @param list_of_temp_files: An iterable of strings that indicate file locations to be deleted.
    @param deleteMe: determines if the files are actually deleted 
    """
    sys.stderr.write("\n%sFiles deleted:\n" % (logstr))
    for f in list_of_temp_files:
        if deleteMe:
            os.remove(f)
        sys.stderr.write("\t%s\n" % (f))
    
    
def main():
    """
    The program's entry point.  This creates the command line interface, parses the options, and dispatches the processes for alignment.
    """
    now = datetime.datetime.now()
    usage = "usage: %prog [options] <reference_genome_file> <output_file> { <reads1> <reads2> | <reads> }\n\nExample: \n %prog ~/ref/human.fa ~/output/human1000.sam ~/reads/human_single_end.fastq"
    version = "%prog " + Constants.version
    description = "This For single end reads, the file location should be the last argument.  For paired end data, there should be two files.  For hairpin data, the two read files should start with the T-enriched one and be followed by the A-enriched one."
    epilog = Constants.creation_string
    p = optparse.OptionParser(usage = usage, version = version, 
                              description = description, epilog = epilog, 
                              formatter = IndentedHelpFormatterWithNL.IndentedHelpFormatterWithNL())
    protocol_group = optparse.OptionGroup(p, "Options for the layout and the construction protocol.")
    protocol_group.add_option('--protocol', '-p', help='This gives the construction protocol that will be used. [default: %default]\n' + 
                 str(Constants.PROTOCOL_DIRECTIONAL) + ': directional\n' + 
                 str(Constants.PROTOCOL_BIDIRECTIONAL) + ': bidirectional\n' + 
                 str(Constants.PROTOCOL_PBAT) + ': PBAT\n' +
                 str(Constants.PROTOCOL_HAIRPIN) + ': hairpin\n' ,
                 #str(Constants.PROTOCOL_ONLYCT) + ': only CT reads against a CT genome\n' +
                 #str(Constants.PROTOCOL_ONLYGA) + ': only GA reads against a GA genome\n' ,
                 default=Constants.PROTOCOL_DIRECTIONAL)
    protocol_group.add_option('--hairpin_paired', '-H', help = 'Treat unrecovered sequences for hairpin data processing as paired end data rather than single end data.', action='store_true', default = False)
    protocol_group.add_option('--insertSizeAvg', '-v', help='Specifies the mean insert size to use when pairing.', default=None)
    protocol_group.add_option('--insertSizeStdDev', '-d', help='Specifies the standard deviation of the insert size to use when pairing', default=None)
    p.add_option_group(protocol_group)
    file_group = optparse.OptionGroup(p, "File management options.")
    file_group.add_option('--path', '-P', help='The path to the BFAST (or BFAST-Gap) executable file.  If this option is not used, then the PATH variable will be searched for the executable.', default = None)
    file_group.add_option('--bfast_gap', '-g', help='Use the BFAST-Gap executable in the system PATH variable.', action='store_true', default=False)
    file_group.add_option('--tmpDir','-T',help='Specifies the directory to store temporary files. [default: the directory where the outputfile is located.]', default = None)
    file_group.add_option('--keep_sam', '-k', help='When this flag is set, the raw SAM files created by BFAST are kept in the temporary directory. [default: %default]', action='store_true', default=False)
    file_group.add_option('--gzip', '-z', help='The input file is gzip compressed, and the output file will be gzip compressed. [default: %default]', action='store_true', default=False)
    file_group.add_option('--startReadNum', '-s', help='Specifies the read to begin with.  Reads are numbered starting at 0. (skip the first startReadNum-1 reads) [default: %default]', default = None)
    file_group.add_option('--endReadNum', '-e', help='Specifies the last read to use (inclusive) [default: %default]', default = None)
    file_group.add_option('--removeComments', '-C', help='Do not include comments at the top of the SAM file.  This can be used while partitioning the input so that the output SAM files can be easily concatenated.', action='store_true', default=False)
    p.add_option_group(file_group)
    timing_group = optparse.OptionGroup(p, "Thread, process, and time management options.")
    timing_group.add_option('--numAlignmentProcesses', '-j', help = "The number of alignment processes to run at the same time.  This can be used for memory or processing core constrained systems to run fewer alignment processes at once. [default: %default]", default = '4')
    timing_group.add_option('--usePipe', '-u', help = 'Pipe the output of BFAST match into BFAST local_align.  Not piping output uses more secondary memory, but piping output probably uses more primary memory.  This prevents the intermediate and temporary matches file from being created. [default: %default]', action='store_true', default=False)
    timing_group.add_option('--numThreads', '-n', help='The number of threads for each BFAST process to use.  Please be aware that runtime efficiency is adversely impacted if the total number of threads is larger than the available processing cores. [default: %default]', default='1')
    timing_group.add_option('--multiProcessPost', '-W', help='Use multiple processes for post processing.  This multiprocessing uses six processes, so this switch should only be engaged if the machine has six or more cores. [default: %default]', action = 'store_true', default = False)
    p.add_option_group(timing_group)
    alignment_group = optparse.OptionGroup(p, "BFAST matching and alignment options")
    alignment_group.add_option('--keySize', '', help = 'Specifies to truncate all indexes to have the given key size (must be greater than the hash width) [default: %default]', default = None)
    alignment_group.add_option('--maxKeyMatches', '', help = 'Specifies the maximum number of matches to allow before a key is ignored [default: %default]', default = None)
    alignment_group.add_option('--maxNumReadMatches', '', help = 'Specifies the maximum total number of matches to consider before the read is discarded [default: %default]', default = None)
    alignment_group.add_option('--maxNumAlignmentMatches', '', help = 'Specifies the maximum number of candidates to initiate alignment for a given match [default: %default]', default = None)
    alignment_group.add_option('--mainIndexes', '-i', help = 'The index numbers for the main bif files (comma separated) [default: %default]', default = None)
    alignment_group.add_option('--secondaryIndexes', '-I', help = 'The index numbers for the secondary bif files (comma separated) [default: %default]', default = None)   
    alignment_group.add_option('--offsets', '-o', help = 'Specifies the offsets for matching. [default: %default]', default = None)
    p.add_option_group(alignment_group)
    alignment_score_group = optparse.OptionGroup(p, "BFAST Alignment Score Options","Options for determining the function to use for calculating the alignment score.")
    alignment_score_group.add_option('--gap_open', '', help='The gap open penalty for scoring an alignment. [default: %default]', default = Constants.GAP_OPEN)
    alignment_score_group.add_option('--gap_extension', '', help='The gap extension penalty for scoring an alignment. [default: %default]', default = Constants.GAP_EXTENSION)
    alignment_score_group.add_option('--match_score', '', help='The score value for a matching base for scoring an alignment. [default: %default]', default = Constants.MATCH_SCORE)
    alignment_score_group.add_option('--mismatch_score', '', help='The base mismatch penalty for scoring an alignment. [default: %default]', default = Constants.MISMATCH_SCORE)
    alignment_score_group.add_option('--scoringMatrixFileName', '-x', help='The alignment scoring function can be specified in this file in the format that BFAST uses.  If this is not specified then the alignment score function specified in the command line will be used. [default: %default]', default = None)
    p.add_option_group(alignment_score_group)
    filter_group = optparse.OptionGroup(p, "Filter Options", "Specify the function for filtering out low quality alignments and reporting the reads as unmapped.")
    filter_group.add_option('--cutoff_bs', '-b', help='Function for filtering low quality alignments for bisulfite treated reads. [default: %default]', default=Constants.BISPIN_DEFAULT_CUTOFF_BS)
    filter_group.add_option('--cutoff_r', '-r', help='Function for filtering low quality alignments for recovered reads.  This is for hairpin data only. [default: %default]', default=Constants.BISPIN_DEFAULT_CUTOFF_R)
    p.add_option_group(filter_group)
    rescore_group = optparse.OptionGroup(p, "Rescore Options", "Ambiguously aligned reads can be rescored to see if one will be unique.")
    rescore_group.add_option('--noRescore', '-N', help='Do NOT rescore ambiguously aligned reads to find unique alignments. [default: %default]', action='store_true', default = False)
    rescore_group.add_option('--rescore_matrix_1', '-m', help='The location of the rescoring matrix used to rescore ambiguously aligned reads for a unique score.  If the second rescoring matrix is specified, it is assumed that this rescoring matrix will be for C to T methylation rescoring. The default is:\n' + Constants.RESCORE_DEFAULT, default=None)
    rescore_group.add_option('--rescore_matrix_2', '-M', help='The location of the rescoring matrix used to rescore ambiguously aligned reads for G to A methylation rescoring for a unique score.  [default: %default]', default=None)
    p.add_option_group(rescore_group)
    options, args = p.parse_args()
    if len(args) == 0:
        p.print_help()
    if len(args) < 3 or len(args) > 4:
        p.error('Not enough arguments.  Check the usage.')
    try:
        input_protocol = int(options.protocol)
    except ValueError:
        p.error('The protocol is not an integer.  Please check.')
    if input_protocol < Constants.PROTOCOL_DIRECTIONAL or input_protocol > Constants.PROTOCOL_ONLYGA:
        p.error('The protocol number is invalid.  Please check.')
    rescore = not options.noRescore
    rescore_matrix = BisPin_postprocess.handleRescoreOptions(rescore, options.rescore_matrix_1, options.rescore_matrix_2)
    fastafile = args[0]
    if not os.path.exists(fastafile):
        p.error("The reference genome file does not exist or could not be accessed.")
    outputfile = args[1]
    outdir = os.path.dirname(outputfile)
    try: 
        os.makedirs(outdir)
    except OSError:
        if not os.path.isdir(outdir):
            p.error("The folder for the output could not be created or was unavailable.")
    if not os.path.isdir(outdir):
        p.error("The folder for the output is invalid.  Please check the input string.")
    filename = os.path.basename(fastafile)
    directory = os.path.dirname(fastafile)
    path_to_bfast = None
    if options.path == None and not options.bfast_gap:
        path_to_bfast = BisPin_util.which("bfast")
    elif options.path == None and options.bfast_gap:
        path_to_bfast = BisPin_util.which("bfast-gap")
    if path_to_bfast == None:
        path_to_bfast = BisPin_util.which(options.path)
    if path_to_bfast == None:
        p.error("The BFAST executable could not be found.  Please check the path.")
    numThreads = options.numThreads
    tmpDir = options.tmpDir
    noRC = False #options.noRC
    reads1 = args[2]
    layout = Constants.LAYOUT_SINGLE
    if not os.path.exists(reads1):
        p.error("The reads file could not be found.")
    reads2 = None
    if len(args) == 4:
        reads2 = args[3]
        layout = Constants.LAYOUT_PAIRED
        if not os.path.exists(reads2):
            p.error("The second reads file could not be found.")
    keep_sam = options.keep_sam
    if tmpDir == None:
        tmpDir = os.path.dirname(outputfile)
    if not tmpDir.endswith("/"):
        tmpDir += "/"
    tmpDir = os.path.join(tmpDir)
    if not os.path.isdir(tmpDir):
        p.error("The temporary directory could not be found.")
    try:
        cutoff_bs = float(options.cutoff_bs)
        cutoff_r = float(options.cutoff_r)
    except ValueError:
        p.error("There was a problem in creating the filter cutoff values.")
    if (options.scoringMatrixFileName == None):
        scoring_file = createAlignmentScoringFunctionFile(tmpDir,options.gap_open,options.gap_extension,options.match_score,options.mismatch_score)
    else:
        scoring_file = os.path.join(options.scoringMatrixFileName)
        if not os.path.exists(scoring_file):
            p.error("The specified scoring file %s does not exist.  Please check." % str(scoring_file))
    if (input_protocol == Constants.PROTOCOL_HAIRPIN or layout == Constants.LAYOUT_PAIRED) and (reads2 == None):
        p.error('The input protocol given needs two reads files, but only one was given.')
    start = None if options.startReadNum == None else int(options.startReadNum)
    end = None if options.endReadNum == None else int(options.endReadNum)
    useSecondaryIndexes = False
    if options.secondaryIndexes != None and options.mainIndexes != None:
        useSecondaryIndexes = True
    if (start != None or end != None) and input_protocol == Constants.PROTOCOL_HAIRPIN and useSecondaryIndexes == True:
        p.error("The partitioning feature for the hairpin protocol is not supported at this time for secondary indexes.  The Linux head and tail commands can be used to partition data.")
    if start != None and end != None and (start > end or start < 0 or end < 0):
        p.error("The start read number is either larger than the end read number or they are negative.  Please check the start read number %s and the end read number %s." % (str(start), str(end)))
    manager = multiprocessing.Manager()
    try:
        numAlignmentProcesses = int(options.numAlignmentProcesses)
    except ValueError:
        p.error("The numAlignmentProcesses needs to be an integer.")
    sema = multiprocessing.Semaphore(numAlignmentProcesses)
    usePipe = options.usePipe
    list_of_temp_files = []
    additional_parameters = manager.dict()
    additional_parameters[Constants.PROTOCOL] = str(input_protocol)
    gzip_switch = options.gzip
    removeComments = options.removeComments
    try:
        keySize = int(options.keySize)
    except TypeError:
        keySize = None
    except ValueError:
        p.error("The keySize was not an integer.")
    try:
        maxKeyMatches = int(options.maxKeyMatches)
    except TypeError:
        maxKeyMatches = None
    except ValueError:
        p.error("The maxKeyMatches was not an integer.")
    try:
        maxNumReadMatches = int(options.maxNumReadMatches)
    except TypeError:
        maxNumReadMatches = None
    except ValueError:
        p.error("The maxNumReadMatches was not an integer.")
    try:
        maxNumAlignmentMatches = int(options.maxNumAlignmentMatches)
    except TypeError:
        maxNumAlignmentMatches = None
    except ValueError:
        p.error("The maxNumAlignmentMatches was not an integer.")
    additional_parameters[Constants.GZIP] = gzip_switch
    additional_parameters[Constants.keySize] = keySize
    additional_parameters[Constants.maxKeyMatches] = maxKeyMatches
    additional_parameters[Constants.maxNumReadMatches] = maxNumReadMatches
    additional_parameters[Constants.maxNumAlignmentMatches] = maxNumAlignmentMatches
    additional_parameters[Constants.insertSizeAvg] = options.insertSizeAvg
    additional_parameters[Constants.insertSizeStdDev] = options.insertSizeStdDev
    additional_parameters[Constants.mainIndexes] = options.mainIndexes
    additional_parameters[Constants.secondaryIndexes] = options.secondaryIndexes
    additional_parameters[Constants.offsets] = options.offsets
    additional_parameters[Constants.offsets] = options.offsets
    additional_parameters[Constants.HAIRPIN_PAIRED] = None
    pool = multiprocessing.Pool(processes = Constants.DEFAULTNUMPROCESSES, initializer = init_child, initargs = (sema,))
    processes_get = []
    post_filename = None
    sys.stderr.write("%sStarting BisPin_align.  Current time: %s\n" % (logstr, str(now)))
    command_line_string = "options: %s args: %s rescore_matrix: %s" % (str(options), str(args), str(rescore_matrix))
    sys.stderr.write("%sCommand line -- %s\n" % (logstr, command_line_string))
    sys.stderr.flush()
    if input_protocol == Constants.PROTOCOL_HAIRPIN:  #Process hairpin data
        sys.stderr.write("%sFinding recovered sequences from the hairpin files located at %s and %s ....\n" % (logstr, reads1, reads2))
        sys.stderr.flush()
        convert_now = datetime.datetime.now()
        hairpin_paired = options.hairpin_paired
        additional_parameters[Constants.HAIRPIN_PAIRED] = hairpin_paired
        hairpin_recovery_files = recover(os.path.join(reads1), os.path.join(reads2), tmpDir, gzip_switch, start, end)
        (R_filename, CT_filename, GA_filename, post_filename) = hairpin_recovery_files
        if keep_sam:
            list_of_temp_files += hairpin_recovery_files[1:3]
        else:
            list_of_temp_files += hairpin_recovery_files[0:3]
        sys.stderr.write("%sCreating C-to-T and G-to-A versions of the unrecovered reads ... \n" % (logstr))
        sys.stderr.flush()
        p_CT = pool.apply_async(BisPin_util.convert_seqs, (tmpDir, CT_filename, True, tmpDir, Constants.FASTQ, start, end, gzip_switch))
        p_GA = pool.apply_async(BisPin_util.convert_seqs, (tmpDir, GA_filename, False, tmpDir, Constants.FASTQ, start, end, gzip_switch))
        CT_filename = p_CT.get()[0]
        GA_filename = p_GA.get()[0]
        list_of_temp_files.append(CT_filename)
        list_of_temp_files.append(GA_filename)
        if hairpin_paired:
            reads_file1 = os.path.basename(reads1)
            CT_GA_file = os.path.join(tmpDir, reads_file1 + ".paired.CT_GA")
            makeOneFilePairedEnd(CT_filename, GA_filename, CT_GA_file, gzip_switch, noRC)
            list_of_temp_files.append(CT_GA_file)
        convert_later = datetime.datetime.now()
        sys.stderr.write("%sFinished converting reads.  Now aligning the reads with BFAST...\n" % (logstr))
        sys.stderr.flush()
        align_now = datetime.datetime.now()
        p_R = pool.apply_async(run_align, (path_to_bfast, R_filename, filename, directory, tmpDir, numThreads, Constants.CONV_R_R, scoring_file, layout, additional_parameters, usePipe, True))
        processes_get.append(p_R)
        if hairpin_paired:
            p_CT_GA_CT = pool.apply_async(run_align, (path_to_bfast, CT_GA_file, filename, directory, tmpDir, numThreads, Constants.CONV_CT_GA_CT, scoring_file, layout, additional_parameters, usePipe))
            p_CT_GA_GA = pool.apply_async(run_align, (path_to_bfast, CT_GA_file, filename, directory, tmpDir, numThreads, Constants.CONV_CT_GA_GA, scoring_file, layout, additional_parameters, usePipe))
            processes_get.append(p_CT_GA_CT)
            processes_get.append(p_CT_GA_GA)
        else:
            p_CT_CT = pool.apply_async(run_align, (path_to_bfast, CT_filename, filename, directory, tmpDir, numThreads, Constants.CONV_CT_CT, scoring_file, layout, additional_parameters, usePipe))
            p_GA_GA = pool.apply_async(run_align, (path_to_bfast, GA_filename, filename, directory, tmpDir, numThreads, Constants.CONV_GA_GA, scoring_file, layout, additional_parameters, usePipe))
            processes_get.append(p_CT_CT)
            processes_get.append(p_GA_GA)
    elif layout == Constants.LAYOUT_PAIRED: #Process paired end data #directional -- CTread1GAread2CTgenome and CTread1GAread2GAgenome
        sys.stderr.write("%sCreating converted versions of the reads located at %s and %s ...\n" % (logstr, reads1, reads2))
        sys.stderr.flush()
        convert_now = datetime.datetime.now()
        reads_file1 = os.path.basename(reads1)
        reads_dir1 = os.path.dirname(reads1)
        reads_file2 = os.path.basename(reads2)
        reads_dir2 = os.path.dirname(reads2)
        CT_GA_file = os.path.join(tmpDir, reads_file1 + ".paired.CT_GA")
        GA_CT_file = os.path.join(tmpDir, reads_file1 + ".paired.GA_CT")
        if input_protocol == Constants.PROTOCOL_PBAT:
            p_GA = pool.apply_async(BisPin_util.convert_seqs, (reads_dir1, reads_file1, False, tmpDir, Constants.FASTQ, start, end, gzip_switch))
            p_CT = pool.apply_async(BisPin_util.convert_seqs, (reads_dir2, reads_file2, True, tmpDir, Constants.FASTQ, start, end, gzip_switch))
            paired1 = p_GA.get()[0]
            paired2 = p_CT.get()[0]
            makeOneFilePairedEnd(paired1, paired2, GA_CT_file, gzip_switch, noRC)
            list_of_temp_files += [paired1, paired2]
            list_of_temp_files.append(GA_CT_file)
        elif input_protocol == Constants.PROTOCOL_BIDIRECTIONAL:
            p_GA1 = pool.apply_async(BisPin_util.convert_seqs, (reads_dir1, reads_file1, False, tmpDir, Constants.FASTQ, start, end, gzip_switch))
            p_CT1 = pool.apply_async(BisPin_util.convert_seqs, (reads_dir1, reads_file1, True, tmpDir, Constants.FASTQ, start, end, gzip_switch))
            p_GA2 = pool.apply_async(BisPin_util.convert_seqs, (reads_dir2, reads_file2, False, tmpDir, Constants.FASTQ, start, end, gzip_switch))
            p_CT2 = pool.apply_async(BisPin_util.convert_seqs, (reads_dir2, reads_file2, True, tmpDir, Constants.FASTQ, start, end, gzip_switch))
            paired1 = p_CT1.get()[0]
            paired2 = p_GA2.get()[0]
            paired3 = p_GA1.get()[0]
            paired4 = p_CT2.get()[0]
            p_1 = pool.apply_async(makeOneFilePairedEnd, (paired1, paired2, CT_GA_file, gzip_switch, noRC))
            p_2 = pool.apply_async(makeOneFilePairedEnd, (paired3, paired4, GA_CT_file, gzip_switch, noRC))
            p_1.get()
            p_2.get()
            list_of_temp_files += [paired1, paired2, paired3, paired4]
            list_of_temp_files.append(CT_GA_file)
            list_of_temp_files.append(GA_CT_file)
        else:
            p_CT = pool.apply_async(BisPin_util.convert_seqs, (reads_dir1, reads_file1, True, tmpDir, Constants.FASTQ, start, end, gzip_switch))
            p_GA = pool.apply_async(BisPin_util.convert_seqs, (reads_dir2, reads_file2, False, tmpDir, Constants.FASTQ, start, end, gzip_switch))
            paired1 = p_CT.get()[0]
            paired2 = p_GA.get()[0]
            makeOneFilePairedEnd(paired1, paired2, CT_GA_file, gzip_switch, noRC)
            list_of_temp_files += [paired1, paired2]
            list_of_temp_files.append(CT_GA_file)
        convert_later = datetime.datetime.now()
        sys.stderr.write(logstr + "Finished converting reads.  Now aligning reads with BFAST...\n")
        sys.stderr.flush()
        align_now = datetime.datetime.now()
        conversions = ""
        if input_protocol == Constants.PROTOCOL_DIRECTIONAL or input_protocol == Constants.PROTOCOL_BIDIRECTIONAL: # CT_GA_GA directional
            conversions += "CT_GA_GA "
            p_CT_GA_GA = pool.apply_async(run_align, (path_to_bfast, CT_GA_file, filename, directory, tmpDir, numThreads, Constants.CONV_CT_GA_GA, scoring_file, layout, additional_parameters, usePipe))
            processes_get.append(p_CT_GA_GA)
        if input_protocol == Constants.PROTOCOL_DIRECTIONAL or input_protocol == Constants.PROTOCOL_BIDIRECTIONAL: #CT_GA_CT directional
            conversions +=  "CT_GA_CT "
            p_CT_GA_CT = pool.apply_async(run_align, (path_to_bfast, CT_GA_file, filename, directory, tmpDir, numThreads, Constants.CONV_CT_GA_CT, scoring_file, layout, additional_parameters, usePipe))
            processes_get.append(p_CT_GA_CT)
        if input_protocol == Constants.PROTOCOL_BIDIRECTIONAL or input_protocol == Constants.PROTOCOL_PBAT: #GA_CT_CT pbat
            conversions += "GA_CT_CT "
            p_GA_CT_CT = pool.apply_async(run_align, (path_to_bfast, GA_CT_file, filename, directory, tmpDir, numThreads, Constants.CONV_GA_CT_CT, scoring_file, layout, additional_parameters, usePipe))
            processes_get.append(p_GA_CT_CT)
        if input_protocol == Constants.PROTOCOL_BIDIRECTIONAL or input_protocol == Constants.PROTOCOL_PBAT: # GA_CT_GA pbat
            conversions += "GA_CT_GA "
            p_GA_CT_GA = pool.apply_async(run_align, (path_to_bfast, GA_CT_file, filename, directory, tmpDir, numThreads, Constants.CONV_GA_CT_GA, scoring_file, layout, additional_parameters, usePipe))
            processes_get.append(p_GA_CT_GA)
        sys.stderr.write(logstr + "Aligning with the following coverted types: " + conversions)
    else: # single end data
        reads_file = os.path.basename(reads1)
        reads_dir = os.path.dirname(reads1)
        sys.stderr.write("%sCreating converted versions of the reads located at %s ...\n" % (logstr, reads1))
        sys.stderr.flush()
        convert_now = datetime.datetime.now()
        if input_protocol == Constants.PROTOCOL_ONLYGA or input_protocol == Constants.PROTOCOL_PBAT:
            p_GA = pool.apply_async(BisPin_util.convert_seqs, (reads_dir, reads_file, False, tmpDir, Constants.FASTQ, start, end, gzip_switch))
            GA_filename = p_GA.get()[0]
            list_of_temp_files.append(GA_filename)
        elif input_protocol == Constants.PROTOCOL_BIDIRECTIONAL:
            p_GA = pool.apply_async(BisPin_util.convert_seqs, (reads_dir, reads_file, False, tmpDir, Constants.FASTQ, start, end, gzip_switch))
            p_CT = pool.apply_async(BisPin_util.convert_seqs, (reads_dir, reads_file, True, tmpDir, Constants.FASTQ, start, end, gzip_switch))
            GA_filename = p_GA.get()[0]
            list_of_temp_files.append(GA_filename)
            CT_filename = p_CT.get()[0]
            list_of_temp_files.append(CT_filename)
        else:
            p_CT = pool.apply_async(BisPin_util.convert_seqs, (reads_dir, reads_file, True, tmpDir, Constants.FASTQ, start, end, gzip_switch))
            CT_filename = p_CT.get()[0]
            list_of_temp_files.append(CT_filename)
        convert_later = datetime.datetime.now()
        sys.stderr.write(logstr + "Finished converting reads.  Now aligning reads with BFAST...\n")
        sys.stderr.flush()
        align_now = datetime.datetime.now()
        if input_protocol == Constants.PROTOCOL_ONLYGA or input_protocol == Constants.PROTOCOL_PBAT or input_protocol == Constants.PROTOCOL_BIDIRECTIONAL: # GA_GA
            p_GA_GA = pool.apply_async(run_align, (path_to_bfast, GA_filename, filename, directory, tmpDir, numThreads, Constants.CONV_GA_GA, scoring_file, layout, additional_parameters, usePipe))
            processes_get.append(p_GA_GA)
        if input_protocol == Constants.PROTOCOL_BIDIRECTIONAL or input_protocol == Constants.PROTOCOL_PBAT: #GA_CT
            p_GA_CT = pool.apply_async(run_align, (path_to_bfast, GA_filename, filename, directory, tmpDir, numThreads, Constants.CONV_GA_CT, scoring_file, layout, additional_parameters, usePipe))
            processes_get.append(p_GA_CT)
        if input_protocol == Constants.PROTOCOL_BIDIRECTIONAL or input_protocol == Constants.PROTOCOL_ONLYCT or input_protocol == Constants.PROTOCOL_DIRECTIONAL: #CT_CT
            p_CT_CT = pool.apply_async(run_align, (path_to_bfast, CT_filename, filename, directory, tmpDir, numThreads, Constants.CONV_CT_CT, scoring_file, layout, additional_parameters, usePipe))
            processes_get.append(p_CT_CT)
        if input_protocol == Constants.PROTOCOL_BIDIRECTIONAL or input_protocol == Constants.PROTOCOL_DIRECTIONAL: # CT_GA
            p_CT_GA = pool.apply_async(run_align, (path_to_bfast, CT_filename, filename, directory, tmpDir, numThreads, Constants.CONV_CT_GA, scoring_file, layout, additional_parameters, usePipe))
            processes_get.append(p_CT_GA)
    sam_files_list = []
    alignment_phase_list = []
    for p in processes_get:
        (process_files, timing_info) = p.get()
        alignment_phase_list.append(timing_info)
        sam_files_list.append(process_files[0])
        if keep_sam:
            process_files = process_files[1:len(process_files)]
        list_of_temp_files += process_files
    if options.scoringMatrixFileName == None:
        list_of_temp_files.append(scoring_file)
    if keep_sam and input_protocol == Constants.PROTOCOL_HAIRPIN:
            sam_files_list.append(post_filename)
    elif not keep_sam and input_protocol == Constants.PROTOCOL_HAIRPIN:
            list_of_temp_files.append(post_filename)
    align_later = datetime.datetime.now()
    sys.stderr.write(logstr + "Finished aligning reads with BFAST.  Now selecting the best alignments and producing the SAM file(s) ...\n")
    sys.stderr.flush()
    pool.terminate()
    pool = None
    post_now = datetime.datetime.now()
    read_files = {
                  Constants.FASTA: fastafile,
                  Constants.BISPIN_CUTOFF_BS: cutoff_bs, 
                  Constants.BISPIN_CUTOFF_R: cutoff_r,
                  Constants.BISPIN_READS1: reads1,
                  Constants.BISPIN_READS2: reads2,
                  Constants.BISPIN_RESCORE_MATRIX: rescore_matrix,
                  Constants.GZIP: gzip_switch,
                  Constants.STARTREADNUM : start,
                  Constants.ENDREADNUM : end,
                  Constants.SINGLEPROCESSPOST : not options.multiProcessPost,
                  Constants.USESECONDARYINDEXES : useSecondaryIndexes,
                  Constants.HAIRPINPOSTFILE : post_filename,
                  Constants.REMOVECOMMENTS : removeComments,
                  }
    report_string, genome_counstruction_time = BisPin_postprocess.performPostProcessFromMain(sam_files_list, 
                                                   pool, 
                                                   manager,
                                                   outputfile, 
                                                   input_protocol, 
                                                   layout,
                                                   command_line_string, 
                                                   read_files)
    sys.stderr.write("%sDeleting temporary files." % (logstr))
    sys.stderr.flush()
    deleteTempFiles(list_of_temp_files)
    post_later = datetime.datetime.now()
    later = datetime.datetime.now()
    timing_report = BisPin_util.createTimingInfoReport(now, later, alignment_phase_list, 
                                                        align_later - align_now, convert_later - convert_now, 
                                                        post_now, post_later, genome_counstruction_time)
    final_report = report_string + "\n" + timing_report + "\n\n" + command_line_string
    sys.stdout.write(final_report + "\n")
    open(os.path.join(outputfile + '.BisPin.report'), 'w').write(final_report)
    sys.stderr.write("%sEnding BisPin_align.  Current time: %s Elapsed time: %s \n" % (logstr, str(later), str(later - now)))
    
    
if __name__ == '__main__':
    main()