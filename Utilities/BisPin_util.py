import os
import math
import datetime
import gzip
from Utilities import SeqIterator
from Utilities import Constants
#from Utilities import Deprecator

"""
@author: Jacob Porter
@summary: A collection of useful functions shared by the alignment and postprocessing modules.
@see: BisPin_align and BisPin_postprocess
"""


def extractPairedID(original_id1, original_id2):
    """
    This function takes two record sequence IDs and returns their overlapping segments.
    This is for use with paired end data where the sequence IDs have a different end part.
    This occurs with Sherman simulated data.
    @param: original_id1: one sequence id
    @param: original_id2: another sequence id
    @return: The overlapping segments of the sequence ids.
    """
    if original_id1 == None or original_id2 == None:
        return None
    original_id1 = original_id1.strip().split()[0]
    original_id2 = original_id2.strip().split()[0]
    if original_id1 == original_id2:
        return original_id1
    else:
        new_id = "".join(a for a, b in zip(original_id1, original_id2) if a == b)
        if new_id == "":
            raise ValueError("The paired end sequences ids did not overlap.  The FASTQ ids for paired end data should be in the same order and should have similar ids: " + original_id1 + " and " + original_id2)
        return new_id
        

def getPossibleGZIPFile(outputfile, gzip_switch, addGZ = False, mode = 'w'):
    """
    Gets a file object that can be either gzipped or not.
    @param outputfile:
    @param gzip_switch:
    @param addGZ:
    @param mode:
    @return A file object for writing to.
    """
    if gzip_switch:
        if not outputfile.endswith(".gz") and addGZ:
            outputfile += ".gz"
        fd = gzip.open(outputfile, mode + 'b')
    else:
        fd = open(outputfile, mode)
    return fd


def reverseComplement(dna_sequence):
    """
    Computes the reverse complement of a DNA sequence
    @param dna_sequence:  The dna sequence to reverse
    @return: a string representing the reverse complement of the dna sequence 
    """
    complements = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(complements.get(base, base) for base in reversed(dna_sequence))
    

def getReadStrFastaStr(type_str):
    """
    Given a conversion type string, 
    a string representing the read conversion and the genome conversion is returned.
    Used by the alignment processes.
    @see:  See Constants for examples.
    @param type_str: A string representing the read and genome conversions
    @return: The read conversion and genome conversion information.
    """
    type_split = type_str.split(Constants.CONV_SEP)
    fasta_str = type_split[len(type_split) - 1]
    if len(type_split) == 3:
        read_str = type_split[0] + Constants.CONV_SEP + type_split[1]
    else:
        read_str = type_split[0]
    return (read_str, fasta_str)


def changeFastqRecordToSAMRecord(fastq_record, seq_id = None):
    """
    Changes a FASTQ record into a SAM record.  The FLAG is set to unmapped.
    Alignment information does not exist, so it will not be included
    @param fastq_record: A FASTQ tuple record to convert.
    @param seq_id: The sequence ID to use for QNAME.
    """
    sam_dictionary = {}
    if seq_id == None:
        sam_dictionary[Constants.SAM_KEY_QNAME] = fastq_record[0]
    else:
        sam_dictionary[Constants.SAM_KEY_QNAME] = seq_id
    sam_dictionary[Constants.SAM_KEY_FLAG] = Constants.SAM_VALUE_UNMAPPED
    sam_dictionary[Constants.SAM_KEY_RNAME] = Constants.SAM_VALUE_STAR
    sam_dictionary[Constants.SAM_KEY_POS] = Constants.SAM_VALUE_ZERO
    sam_dictionary[Constants.SAM_KEY_MAPQ] = Constants.SAM_VALUE_ZERO
    sam_dictionary[Constants.SAM_KEY_CIGAR] = Constants.SAM_VALUE_STAR
    sam_dictionary[Constants.SAM_KEY_RNEXT] = Constants.SAM_VALUE_STAR
    sam_dictionary[Constants.SAM_KEY_PNEXT] = Constants.SAM_VALUE_ZERO
    sam_dictionary[Constants.SAM_KEY_TLEN] = Constants.SAM_VALUE_ZERO
    sam_dictionary[Constants.SAM_KEY_SEQ] = fastq_record[1]
    sam_dictionary[Constants.SAM_KEY_QUAL] = fastq_record[2]
    sam_dictionary[Constants.SAM_KEY_ALIGNMENT_SCORE] = Constants.SAM_VALUE_BFAST_UNMAPPED_SCORE
    sam_dictionary[Constants.SAM_KEY_PROGRAM] = Constants.SAM_VALUE_PROGRAM
    return sam_dictionary


def tokenizeCigar(cigar):
    """
    This function tokenizes the cigar string for easier processing.
    For example, the string '2M1I27M1D11M' is returned as the list of tuples:
    [(2, 'M'), (1, 'I'), (27, 'M'), (1, 'D'), (11, 'M')]
    @param cigar: The CIGAR string from an alignment.
    @return: A list representing the tokenized CIGAR string.
    """
    cigar = cigar.replace(" ", "")
    tokenized_list = []
    a_number = ""
    for c in cigar:
        if c.isdigit():
            a_number += c
        else:
            tokenized_list.append((int(a_number), c))
            a_number = ""
    return tokenized_list
            
            
def tokenizeMDtag(mdtag):
    """
    This function tokenizes the MD tag string.
    If the string is '5T1T0T21G1^A2' the list is returned:
    [5, 'T', 1, 'T', 0, 'T', 21, 'G', 1, '^A', 2]
    @param mdtag: The MD tag string
    @return: A list representing the MD tag.
    """
    mdtag = mdtag.replace(" ", "")
    tokenized_list = []
    a_number = ""
    a_letter = ""
    if len(mdtag) == 0:
        return tokenized_list
    atdigit = mdtag[0].isdigit()
    for c in mdtag:
        if c.isdigit():
            a_number += c
        else:
            a_letter += c
        if atdigit and not c.isdigit():
            tokenized_list.append(int(a_number))
            atdigit = False
            a_number = ""
        elif not atdigit and c.isdigit():
            tokenized_list.append(a_letter)
            atdigit = True
            a_letter = ""
    if a_number == "":
        tokenized_list.append(a_letter)
    else:
        tokenized_list.append(int(a_number))
    return tokenized_list


def findDeletions(cigar_token, md_token, reference_sequence):
    """
    Finds all the deletions in the alignment.
    @param cigar_token: The tokenized CIGAR string.
    @param md_token: The tokenized MD tag string.
    @param reference_sequence: 
    @return: A list of all the deletions.
    """
    current_position = 0
    deletions_length = 0
    deletion_list = []
    for token in md_token:
        if isinstance(token, int):
            current_position += token
        elif token.startswith('^'):
            reference_index = current_position + deletions_length
            new_delete = reference_sequence[reference_index : reference_index + len(token) - 1]
            deletion_list.append((current_position+1, '^' + new_delete))
            deletions_length += len(token) - 1
        elif isinstance(token, basestring):
            current_position += len(token)
    return deletion_list
            

def averagePhredStrings(phred1, phred2, offset = 33):
    """
    Averages the two phred strings phred1 and phred2 encoded with offset.
    This is used by the hairpin recovery process and gives the PHRED string
    for hairpin recovered strings. 
    @param phred1: A PHRED quality string.
    @param phred2: A PHRED quality string.
    @param offset: The ASCII offset that encodes the PHRED strings.  
    @return: A PHRED quality string that represents the average of two other PHRED strings. 
    """
    return "".join([averagePhredCharacter(c[0], c[1], offset=offset) for c in zip(phred1, phred2) ])


def averagePhredCharacter(c1, c2, offset = 33):
    """
    Average the two phred characters.
    Formula for converting Q score into probability (phred 33): 10**((ord('A')-33)/(-10.0))
    Q = -10*math.log10(P)
    P = 10**(-Q/10)
    opposite of ord() is chr()
    @param c1: A phred character.
    @param c2: Another phred character.
    @param offset: The offset to use for computing the phred stuff.
    @return: The phred character representing the average of c1 and c2.
    """
    p1 = 10**((ord(c1)-offset)/(-10.0)) #Converts a phred character to a probability
    p2 = 10**((ord(c2)-offset)/(-10.0))
    pa = (p1 + p2)/2.0
    q = int(round(-10*math.log10(pa)))
    if offset==33 and q > 41:
        q = 41
    elif q > 40:
        q = 40
    elif q < 0:
        q = 0
    return chr(q+offset)


def convert_seqs(directory, filename, CTorGA, tmpDir = None, file_type = Constants.FASTA, 
                 start = None, end = None, gzip_switch = False, outputfile = None, reverse = False):
    """
    Does the C to T or the G to A conversion of DNA sequences.
    @param directory: A directory where the file to convert will be found.
    @param filename: The filename of the file to convert.
    @param CTorGA: A boolean value when True does a C to T conversion.  Otherwise, a G to A conversion.
    @param tmpDir: A directory to write the converted file to.
    @param file_type: The file type of the converted records.  i.e., FASTA or FASTQ
    @param start: The record number in the sequence of records to start with.
    @param end: The record number in the sequence of records to end with.
    @param gzip_switch: If True, the output is written in gzipped compressed format.
    The input must be gzipped as well.
    @param outputfile: The location of the outputfile to write to.
    @param reverse: Reverses the newly converted sequence when True.
    @return: A list consisting of the output file represented as a string.
    """
    my_directory = directory if tmpDir == None else tmpDir
    if outputfile == None:
        str1 = Constants.CONV_CT if CTorGA else Constants.CONV_GA
        outputfile = os.path.join(my_directory, filename + ".BisPin." + str1)
    fd_convert = getPossibleGZIPFile(outputfile, gzip_switch)
    unconverted_seqs = SeqIterator.SeqIterator(os.path.join(directory, filename), 
                                               file_type=file_type, gzip_switch = gzip_switch)
    converted_seqs = SeqIterator.SeqWriter(fd_convert, file_type=file_type)
    counter = -1
    for rec in unconverted_seqs:
        counter += 1
        if (start != None and counter < start) or (end != None and counter > end):
            continue
        my_seq = rec[1]
        new_seq = ""
        for base in my_seq:
            if CTorGA:
                if base == "C":
                    base = "T"
                elif base == "c":
                    base = "t"
            else:
                if base == "G":
                    base = "A"
                elif base == "g":
                    base = "a"
            new_seq += base
        if reverse:
            new_seq = new_seq[::-1]
        if file_type == Constants.FASTQ:
            my_seq = (rec[0], new_seq, rec[2], rec[3])
        else:
            my_seq = (rec[0], new_seq)
        converted_seqs.write(my_seq)
    fd_convert.close()
    return [outputfile]


def createTimingInfoReport(now, later, alignment_phase_list, 
                           alignment_time_total, conversion_time, 
                           postprocessing_now, postprocessing_later, genome_counstruction_time):
    """
    Creates the timing information report.
    The parameters are meant to be datetime.timedelta objects.
    @param total_time: The total elapsed wall time from beginning to end
    @param alignment_phase_list: A list of dictionaries of timedeltas representing each  
    phase of each alignment process
    @param alignment_time_total: The total wall time that alignment took
    @param conversion_time: The time it took to convert the reads
    @param postprocessing_time: The wall time it took for BisPin to do its post processing 
    @param genome_counstruction_time: The amount of time BisPin took to load the reference genome 
    into memory for its post processing.
    @return: A string representing the timing information report.
    """
    phase_time = {}
    counter = {}
    if alignment_phase_list != None:
        for d in alignment_phase_list:
            for key in d:
                phase_time[key] = phase_time.get(key, datetime.timedelta(0)) + d[key]
                counter[key] = counter.get(key, 0) + 1
    for key in phase_time:
        phase_time[key] = phase_time[key] / counter[key]
    timing_report = ""
    timing_report += "Timing Report\n"
    timing_report += "-------------\n"
    if now != None and later != None:
        total_time = later - now
        timing_report += "The process was started at %s and ended at %s.\n" % (str(now), str(later))
    if postprocessing_now != None and postprocessing_later != None:
        postprocessing_time = postprocessing_later - postprocessing_now
        timing_report += "The BisPin post processing was started at %s and ended at %s.\n" % (str(postprocessing_now), str(postprocessing_later))
    if now != None and later != None:
        timing_report += "\nTime for the whole process:\t%s\n" % str(total_time)
    if conversion_time != None:
        timing_report += "Time for converting the reads:\t%s\n" % str(conversion_time)
    if alignment_time_total != None:
        timing_report += "Time for the alignment phase:\t%s" % str(alignment_time_total)
    if postprocessing_now != None and postprocessing_later != None:
        timing_report += "\nTime for BisPin post processing:\t%s\n" % str(postprocessing_time)
        if genome_counstruction_time != None:
            timing_report += "\tTime to load the reference genome:\t%s\n" % str(genome_counstruction_time)
            timing_report += "\tTime to modify the SAM records:\t%s\n" % str(postprocessing_time - genome_counstruction_time)
    if alignment_phase_list != None:
        timing_report += "\n"
        timing_report += "Alignment phase average timing breakdown:\n"
        for key in phase_time:
            timing_report += "\t%s\t%s\n" % (str(key), str(phase_time[key]))
    return timing_report


def which(program):
    """
    Code taken from:
    http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python/377028#377028
    This checks if an executable is present.
    @param program: Either a path to a program or the name of a program on the system PATH
    @return: The path to the program.  If the program does not exist, return None.
    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, _ = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

