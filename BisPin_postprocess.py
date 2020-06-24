#!/usr/bin/python
import sys
from functools import reduce
sys.dont_write_bytecode = True
import optparse
import datetime
import os
import multiprocessing
from Utilities import SeqIterator
from Utilities import SeqDoubleIterator
from Utilities import Constants
from Utilities import BisPin_util
from Utilities import IndentedHelpFormatterWithNL
#from Utilities import Deprecator

"""
@author: Jacob Porter
@summary: The postprocessing module selects the best scoring alignments from BFAST SAM files and modifies the flag and the MD tag as appropriate.
@requires:  The SAM files (from BisPin_align) must exist and must use the correct naming convention. 
TODO: Make sure that the postprocessing can handle reference and reads with wildcard N characters.
TODO: BFAST gives incorrect flags for paired end mapping that indicate incorrect orientations.
            #These flags are for the GA converted genome
            #115 should be 83 RC
            #179 should be 147 not RC
            #These flags are for the CT converted genome
            #65 should be 99 not RC
            #129 should be 147 RC
TODO: Bismark has a bug in PBAT alignments.  When the reverse complement is mapped, the SAM flag does not indicate this.
TODO: The program barfs on the multiple index functionality when using a partition and multiprocessing and multiple indexes.  This appears to be fixed.
TODO: Test the new rescoring functionality for positive C to T matrix scores.
TODO: Does the rescoring correctly score for the right direction when it looks at methylation calls?
TODO: Hairpin recovered methylation calling statistics are substantially different compared to not using recovery.  This could be a bug.  
TODO: Using BFAST-Gap there is sometimes a read that cannot be found in the alignmments file.  This could be because of multiple indexes. Currently, BisPin skips past records that it cannot find.
TODO: Perhaps the secondary index functionality should always be used in case there are missing or out of order records
TODO: When a deletion occurs before an insertion, the length of the deletion is shorter in the CIGAR string compared to the MD tag string.  
This occurs rarely with BFAST-Gap with the logistic and exponential settings.  It is unclear if this is related to the changes in BFAST-Gap or if it is a rare bug in BFAST.  The bug is not always present.
Looking at some examples, it makes more sense to prefer a shorter deletion (by one) since this will result in fewer downstream mismatches and a higher score. The CIGAR string is correct, but the MD tag is not.
"""

logstr = "\nBisPin_postprocess: "

def writeSamComments(sam_output, args):
    """
    Writes comments into a SAM file.
    @param sam_output:  A file object to write to.
    @param args: A list/tuple of strings.  The list/tuple can contain lists or tuples of strings.  These represent the comments.
    """
    for arg in args:
        if isinstance(arg, list) or isinstance(arg, tuple):
            for s in arg:
                sam_output.write(str(s))
        else:
            sam_output.write(str(arg))
            

def selectBestAlignmentAndPostProcessSAMRecord(sam_files_list, pool, input_protocol, layout, read_iterator, 
                                               outputfile, cutoff_bs, cutoff_r, rescore_matrix, start, end, 
                                               single_process_switch, gzip_switch, command_line_string, 
                                               manager, genome_file, use_secondary_indexes, hairpin_post, remove_comments):
    """
    This function selects and processes the best alignments from the BFAST alignments.  
    A SAM file is written that indicates the best alignments.
    @param sam_files_list: A list of raw BFAST SAM files to post process.
    @param pool: A pool object for multiprocessing. (Not used.  Set to None.)
    @param input_protocol: A number from Constants indicating the input protocol.
    @param layout: A number from Constants indicating the type of layout.
    @param read_iterator: A single (single end data) or double iterator of the input FASTQ records.  See the SeqIterator file. 
    @param outputfile: A string of the output file for writing the post processed SAM records.
    @param cutoff_bs: The minimum score threshold.  Reads not meeting this score will be marked as filtered.
    @param cutoff_r: The minimum score threshold for recovered reads (from hairpin data).
    @param rescore_matrix: A 2-tuple of dictionaries representing the rescoring function for multireads.  If both elements are set to None, no rescoring will be done.
    @param start: An integer of the first read to process from the input.  If start and end are set to None, process all input.
    @param end: An integer of the last read to process from the input.
    @param single_process_switch: A boolean when set to True will cause only a single process to be used for post-processing.
    @param gzip_switch: A boolean indicating whether the output should be gzip compressed.  If True, compression is turned on.
    @param command_line_string: A string representing the command line parameters and arguments.  Added to the SAM file output.
    @param manager: A multiprocessing manager server for coordinating multiprocessing 
    @param genome_file: A string indicating the location of the genome file. 
    @param use_secondary_indexes: A boolean that is set to True if the SAM files were constructed using secondary indexes.  False otherwise.
    @param hairpin_post: The file location for the recovered reads from hairpin sequencing
    @param remove_comments:  If set to true, no comments will be written to the SAM file  
    @return: A 2-tuple.  In the first half are stats for the alignment efficiency, a breakdown of how reads uniquely mapped, and methylation context counting.  The second half gives the time for loading the reference genome. 
    """
    #Initialize counters and get SAM file information and iterators
    (sam_iterator_dictionary, sam_conversion_info) = makeSamFileTypeDictionary(sam_files_list)
    if use_secondary_indexes:
        (sam_iterator_dictionary2, _) = makeSamFileTypeDictionary(sam_files_list)
    sam_iterator_keys = list(sam_iterator_dictionary.keys())
    filtered_unmapped = 0
    unaligned_unmapped = 0
    ambiguous = 0
    unique = 0
    record_count = 0
    counter = -1
    ambiguous_uniquified = 0
    methylation_counter = {}
    uniquely_mapped_counter = {}
    sys.stderr.write(logstr + "Loading the reference genome into memory.\n")
    sys.stderr.flush()
    #Load the reference genome into memory.
    if single_process_switch:
        (genome_dictionary, id_length_list, genome_construction_time) = makeReferenceGenomeDictionary(genome_file)
        comments = ("@HD\tVN:1.3\tSO:unsorted\tGO:none\n", ["@SQ\tSN:%s\tLN:%s\n" % x for x in id_length_list], 
                    "@PG\tID:%s\tCL:%s\tVN:%s\n" % (Constants.SAM_VALUE_PROGRAM, command_line_string ,Constants.version))
        sam_output = SeqIterator.SeqWriter(BisPin_util.getPossibleGZIPFile(outputfile, gzip_switch, addGZ = True), file_type = Constants.SAM)
        if not remove_comments:
            writeSamComments(sam_output, comments)
    else:
        processes = []
        pipes_parent = []
        pipes_child = []
        return_info_objects = manager.list()
        comments = manager.list()
        q_read = manager.Queue()
        q_write = manager.Queue()
        q_coordinate = manager.Queue()
        for _ in range(Constants.DEFAULTNUMPROCESSES - 1):
            pipe_worker = multiprocessing.Pipe()
            pipes_parent.append(pipe_worker[0])
            pipes_child.append(pipe_worker[1])
        now = datetime.datetime.now()
        finished = manager.Event() #This Event signals when the process is done loading the reference genome.
        #Process for loading the reference genome
        p_genome = multiprocessing.Process(target = multiprocessGenomeService, args = (genome_file, pipes_parent, 
                                                                                       finished, return_info_objects, 
                                                                                       comments, command_line_string))
        p_genome.start()
        finished.wait()
        later = datetime.datetime.now()
        #Start worker processes that will post process the SAM records.
        for i in range(Constants.DEFAULTNUMPROCESSES - 1):
            p_worker = multiprocessing.Process(target = multiprocessPostProcessing, args = (q_read, q_write, 
                                                                                            cutoff_bs, cutoff_r, 
                                                                                            rescore_matrix, input_protocol, 
                                                                                            layout, pipes_child[i]))
            p_worker.start()
            processes.append(p_worker)
        #Start the process that writes the post processed SAM records to the output file.
        p_writer = multiprocessing.Process(target = multiprocessWriteSAMRecord, args = (q_write, outputfile, 
                                                                                        gzip_switch, comments, remove_comments,
                                                                                        Constants.DEFAULTNUMPROCESSES - 1, 
                                                                                        q_coordinate))
        p_writer.start()
        genome_construction_time = later - now
    #Loop through the input reads and post process each one.  This loop assumes that the order of reads
    #by sequence id in the SAM file(s) and the input FASTQ file(s) are the same except when
    #using secondary indexes.  This case is handled with alternate SAM file iterators.
    if input_protocol == Constants.PROTOCOL_HAIRPIN:
        hairpin_post_fd = open(hairpin_post, 'r')
        if use_secondary_indexes:
            hairpin_post_fd_2 = open(hairpin_post, 'r')
    for record in read_iterator:
        counter += 1
        if (start != None and counter < start):
            continue
        if (end != None and counter > end):
            break
        record_count += 1
        write_record = ""
        meilleur_alignement_dictionnaire = {}
        if layout == Constants.LAYOUT_PAIRED or input_protocol == Constants.PROTOCOL_HAIRPIN:
            record1 = record[0]
            record2 = record[1]
            alignement_key = BisPin_util.extractPairedID(record1[0], record2[0])
        else:
            record1 = record
            record2 = None
            record = [record]
            alignement_key = record1[0].split()[0]
        empty = True
        if input_protocol == Constants.PROTOCOL_HAIRPIN:
            hairpin_recovery_location = hairpin_post_fd.readline().split()
            while use_secondary_indexes and hairpin_recovery_location[0] != alignement_key:
                line = hairpin_post_fd_2.readline()
                if line == "" or line == None:
                    break
                hairpin_recovery_location = line.split()
            if hairpin_recovery_location[0] != alignement_key:
                sys.stderr.write(logstr + "The hairpin key from the hairpin recovery post processing file did not match the input FASTQ file: " + hairpin_recovery_location[0] + ", " + alignement_key + "\n")
        #Extract all reads from the SAM files that match the input read(s)
        for key in sam_iterator_keys:
#             if input_protocol == Constants.PROTOCOL_HAIRPIN:
#                 sys.stderr.write("For location %s, we have %s\n" % (str(key), str(hairpin_recovery_location)))
#                 sys.stderr.flush()
            if input_protocol == Constants.PROTOCOL_HAIRPIN and hairpin_recovery_location[1].startswith(Constants.CONV_R_R) and key != Constants.CONV_R_R:
                continue
            elif input_protocol == Constants.PROTOCOL_HAIRPIN and not hairpin_recovery_location[1].startswith(Constants.CONV_R_R) and key == Constants.CONV_R_R:
                continue
            brut_alignements = []
            conversion_info = sam_conversion_info[key]
            next_Id = sam_iterator_dictionary[key].peekAtId()
            if (layout == Constants.LAYOUT_PAIRED or input_protocol == Constants.PROTOCOL_HAIRPIN):
                next_Id = BisPin_util.extractPairedID(alignement_key, next_Id)
#             if input_protocol == Constants.PROTOCOL_HAIRPIN:
#                 sys.stderr.write("Searching for alignement key %s in location %s.  next_Id is: %s\n" % (str(alignement_key), str(key), str(next_Id)))
#                 sys.stderr.flush()
            in_alternate_file = False
            #Check the alternate iterator for the SAM file if secondary indexes are used since 
            #reads that were unmapped in the primary indexes are put at the end of the SAM file.
            if use_secondary_indexes and next_Id != None and not alignement_key == next_Id:
                in_alternate_file = True
                while(True):
                    next_Id = sam_iterator_dictionary2[key].peekAtId()
                    if (layout == Constants.LAYOUT_PAIRED or input_protocol == Constants.PROTOCOL_HAIRPIN):
                        #next_Id = "".join(pairedRegularExpression.split(next_Id))
                        next_Id = BisPin_util.extractPairedID(alignement_key, next_Id)
#                     if input_protocol == Constants.PROTOCOL_HAIRPIN:
#                         sys.stderr.write("The next_Id in the second iterator is: %s\n" % (str(next_Id)))
                    if next_Id != None and alignement_key == next_Id:
#                         if input_protocol == Constants.PROTOCOL_HAIRPIN:
#                             sys.stderr.write("Breaking out of the second iterator.\n")
                        break
                    else:
                        try:
                            next(sam_iterator_dictionary2[key])
                        except StopIteration as si:
                            #If this happens, something is wrong, and the SAM record was not found in the BFAST file.
                            #This could be a bug in BFAST-Gap.
                            sys.stderr.write("The alignment key %s for file type %s could not be found, and the end of the second iterator was reached.\n" % (str(alignement_key), str(key)) )
                            next_Id = None
                            sam_iterator_dictionary2[key].reset()
                            break
                            #raise si
            #while(alignement_key == sam_iterator_dictionary[key].peekAtId()):
            while(next_Id != None and alignement_key == next_Id):
                if in_alternate_file:
                    sam_record = next(sam_iterator_dictionary2[key])
                else:
                    sam_record = next(sam_iterator_dictionary[key])
                if isFirstSegment(sam_record[Constants.SAM_KEY_FLAG]):
                    sam_record[Constants.SAM_KEY_READ] = conversion_info[0]
                    sam_record[Constants.SAM_KEY_GENOME] = conversion_info[2]
                elif isLastSegment(sam_record[Constants.SAM_KEY_FLAG]):
                    sam_record[Constants.SAM_KEY_READ] = conversion_info[1]
                    sam_record[Constants.SAM_KEY_GENOME] = conversion_info[2]
                else:
                    sam_record[Constants.SAM_KEY_READ] = conversion_info[0]
                    sam_record[Constants.SAM_KEY_GENOME] = conversion_info[1]
                brut_alignements.append(sam_record)
                if in_alternate_file:
                    next_Id = sam_iterator_dictionary2[key].peekAtId()
                else:
                    next_Id = sam_iterator_dictionary[key].peekAtId()
                if (layout == Constants.LAYOUT_PAIRED or input_protocol == Constants.PROTOCOL_HAIRPIN):
                    #next_Id = "".join(pairedRegularExpression.split(next_Id))
                    next_Id = BisPin_util.extractPairedID(alignement_key, next_Id)
            if brut_alignements != []:
                empty = False
                meilleur_alignement_dictionnaire[key] = brut_alignements
        #If using only one process, process the read.
        if single_process_switch:
            (write_record, ambiguous_uniquified_integer, direction_id) = computeWriteRecord(meilleur_alignement_dictionnaire, 
                                                                                            empty, 
                                                                                            record1, 
                                                                                            record2, 
                                                                                            alignement_key, 
                                                                                            cutoff_bs, 
                                                                                            cutoff_r, 
                                                                                            rescore_matrix, 
                                                                                            record, 
                                                                                            input_protocol, 
                                                                                            layout, 
                                                                                            methylation_counter, 
                                                                                            genome_dictionary)
            ambiguous_uniquified += ambiguous_uniquified_integer
            if write_record[0] == Constants.SAM_WRITE_UNMAPPED:
                unaligned_unmapped += 1
            elif write_record[0] == Constants.SAM_WRITE_FILTERED:
                filtered_unmapped += 1
            elif write_record[0] == Constants.SAM_WRITE_AMBIGUOUS:
                ambiguous += 1
            elif write_record[0] == Constants.SAM_WRITE_UNIQUE:
                unique += 1
                uniquely_mapped_counter[direction_id] = uniquely_mapped_counter.get(direction_id, 0) + 1
            writeSAMRecord(write_record, sam_output)
        #If using multiple processes, put the reads on a work queue and the read id on another queue
        #so that the output file writer process will write the reads in the same order as the input.
        else:
            q_read.put((meilleur_alignement_dictionnaire, empty, record1, record2, alignement_key, record))
            q_coordinate.put(alignement_key)
    #Once all the work is put onto the work queue, join the processes and extract statistic information.
    if not single_process_switch:
        q_coordinate.put(None)
        for _ in range(Constants.DEFAULTNUMPROCESSES - 1):
            q_read.put(None)
        for proc in processes:
            proc.join()
        p_genome.join()
        for report_info in return_info_objects:
            alignment_distribution = report_info[0]
            unique += alignment_distribution[1]
            ambiguous += alignment_distribution[2]
            unaligned_unmapped += alignment_distribution[3]
            filtered_unmapped += alignment_distribution[4]
            ambiguous_uniquified += alignment_distribution[5]
            p_uniquely_mapped_counter = report_info[1]
            p_methylation_counter = report_info[2]
            for key in p_uniquely_mapped_counter:
                uniquely_mapped_counter[key] = uniquely_mapped_counter.get(key, 0) + p_uniquely_mapped_counter[key]
            for key in p_methylation_counter:
                methylation_counter[key] = methylation_counter.get(key, 0) + p_methylation_counter[key]
        p_writer.join()
    return (([record_count, unique, ambiguous, unaligned_unmapped, filtered_unmapped, ambiguous_uniquified], 
             uniquely_mapped_counter, methylation_counter), genome_construction_time)


def multiprocessGenomeService(genome_file, pipes_parent, finished, return_info_objects, comments, command_line_string):
    """
    For multiprocessing, this function loads the reference genome into a dicitonary and services requests for reference genome info.
    @param genome_file: A string indicating the location of the reference genome FASTA file.
    @param pipes_parent: A list of pipe endpoints.  Worker processes put requests for the reference genome on these pipes.
    @param finished: An event object that indicates when the process is done loading the reference genome into memory.
    @param return_info_objects: A process safe list that is used to report each processes stats.
    @param comments: A multiprocess list for comments.  The SAM writer process will write this to the SAM file.
    @param command_line_string: The command line parameter and argument string to write to the SAM file 
    """
    #Load the reference genome
    (genome_dictionary, id_length_list, _) = makeReferenceGenomeDictionary(genome_file)
    genome_comments = ["@SQ\tSN:%s\tLN:%s\n" % x for x in id_length_list]
    comments.append("@HD\tVN:1.3\tSO:unsorted\tGO:none\n")
    for g in genome_comments:
        comments.append(g)
    comments.append("@PG\tID:%s\tCL:%s\tVN:%s\n" % (Constants.SAM_VALUE_PROGRAM, command_line_string, Constants.version))
    finished.set()
    #Service requests for reference genome information using full duplex pipes.
    while pipes_parent != []:
        remove_pipes = []
        for pipe in pipes_parent:
                    if pipe.poll(0):
                        pipe_object = pipe.recv()
                        if pipe_object[0] == None:
                            remove_pipes.append(pipe)
                            return_info_objects.append(pipe_object[1])
                        else:
                            pipe.send(genome_dictionary[pipe_object[0]][pipe_object[1] : pipe_object[2]])
        for pipe in remove_pipes:
            pipes_parent.remove(pipe)
            pipe.close()
            

def multiprocessPostProcessing(q_read, q_write, cutoff_bs, cutoff_r, rescore_matrix, input_protocol, layout, genome_pipe):
    """
    A function used by multiprocess workers for doing all the post processing of SAM records.  Work is put onto a work queue.
    @param q_read: A multiprocess work queue for extracting work.
    @param q_write: A multiprocess work queue for depositing finished work.
    @param cutoff_bs: A float for filtering low quality reads.
    @param cutoff_r: A float for filtering low quality hairpin recovered reads.
    @param rescore_matrix: A dictionary representing the rescoring function for multireads.  If None, it will not be used.
    @param input_protocol: A number from Constants indicating the input protocol.
    @param layout: A number from Constants indicating the type of layout.
    @param genome_pipe: A pipe to ask for genome reference information.
    """
    #Initializes counters for statistics.
    methylation_counter = {}
    uniquely_mapped_counter = {}
    filtered_unmapped = 0
    unaligned_unmapped = 0
    ambiguous = 0
    unique = 0
    ambiguous_uniquified = 0
    record_count = 0
    #An infinite loop for getting work from a work queue and doing the post processing.
    while(True):
        info = q_read.get()
        if info is None:
            q_write.put(None)
            return_info = ([record_count, unique, ambiguous, unaligned_unmapped, filtered_unmapped, ambiguous_uniquified], 
                           uniquely_mapped_counter, methylation_counter)
            genome_pipe.send((None, return_info))
            return
        else:
            record_count += 1
            (meilleur_alignement_dictionnaire, empty, record1, record2, alignement_key, record) = info
            (write_record, ambiguous_uniquified_integer, direction_id) = computeWriteRecord(meilleur_alignement_dictionnaire, 
                                                                                            empty, 
                                                                                            record1, 
                                                                                            record2, 
                                                                                            alignement_key, 
                                                                                            cutoff_bs, 
                                                                                            cutoff_r, 
                                                                                            rescore_matrix, 
                                                                                            record, 
                                                                                            input_protocol, 
                                                                                            layout, 
                                                                                            methylation_counter, 
                                                                                            None, 
                                                                                            genome_pipe=genome_pipe)
            ambiguous_uniquified += ambiguous_uniquified_integer
            if write_record[0] == Constants.SAM_WRITE_UNMAPPED:
                unaligned_unmapped += 1
            elif write_record[0] == Constants.SAM_WRITE_FILTERED:
                filtered_unmapped += 1
            elif write_record[0] == Constants.SAM_WRITE_AMBIGUOUS:
                ambiguous += 1
            elif write_record[0] == Constants.SAM_WRITE_UNIQUE:
                unique += 1
                uniquely_mapped_counter[direction_id] = uniquely_mapped_counter.get(direction_id, 0) + 1
            q_write.put(write_record)
            
    
def multiprocessWriteSAMRecord(q_write, outputfile, gzip_switch, comments, remove_comments, num_processes, q_coordinate):
    """
    This function writes post processed SAM records to a file with multiprocessing.
    @param q_write: A queue for getting records to write to the SAM file.
    @param outputfile: A string representing the location of the SAM output file.
    @param gzip_switch: A boolean indicating whether the output should be gzip compressed.  If True, compression is turned on.
    @param comments: A multiprocess list for comments.  The SAM writer process will write this to the SAM file.
    @param num_processes: The number of worker processes doing work.  This is used to kill the process when the work is done.
    @param q_coordinate: A coordination data structure that is used to ensure that SAM records have the same order as the input file.
    """
    poison_pills = 0
    write_cache = {} #A cache of SAM records that have been received out of order.
    sam_output = SeqIterator.SeqWriter(BisPin_util.getPossibleGZIPFile(outputfile, gzip_switch, addGZ = True), 
                                       file_type = Constants.SAM)
    if not remove_comments:
        writeSamComments(sam_output, comments)
    sam_output.flush()
    coordinateFull = True
    writerFull = True
    while (coordinateFull and writerFull):
        while (coordinateFull): # An infinite loop that gets the next sequence id to write and checks the write cache if it is there.
            alignement_key = q_coordinate.get()
            if alignement_key == None:
                coordinateFull = False
            else:
                write_cache_record = write_cache.get(alignement_key, None)
                if write_cache_record != None:
                    writeSAMRecord(write_cache_record, sam_output)
                else:
                    break
        while (writerFull): #If a sequence id is not in the write cache, search for it on the queue until it is found.
            write_record = q_write.get()
            if write_record == None:  #Terminates the process if enough 'poison pills' are received. 
                poison_pills += 1
                if poison_pills == num_processes:
                    writerFull = False
            else:
                record_key = write_record[1][0][0][Constants.SAM_KEY_QNAME]
                if record_key.startswith(alignement_key):
                    writeSAMRecord(write_record, sam_output)
                    break
                else:
                    write_cache[record_key] = write_record  #Put records that do not match the record id on the cache.
    sam_output.flush()
    sam_output.close()


def computeWriteRecord(meilleur_alignement_dictionnaire, empty, record1, record2, alignement_key, 
                       cutoff_bs, cutoff_r, rescore_matrix, record, input_protocol, layout, 
                       methylation_counter, genome_dictionary, genome_pipe = None):
    """
    Does the work of post processing a SAM record and generates a write_record object.
    @param meilleur_alignement_dictionnaire: A dictionary representing the SAM records to choose and process.
    @param empty: A boolean indicating if the meilleur_alignement_dictionnaire is empty.
    @param record1: The first pair FASTQ input record for paired end data.  The sole single record for single end data.
    @param record2: The second pair FASTQ record for paired end data (and hairpin data).  None for single end data.
    @param alignement_key: The sequence id of the FASTQ record.  Must match the SAM QNAME field.
    @param cutoff_bs: A float for filtering low quality reads.
    @param cutoff_r: A float for filtering low quality hairpin recovered reads.
    @param rescore_matrix: A 2-tuple of dictionaries representing the rescoring function for multireads.  If both elements are set to None, no rescoring will be done.
    @param record: The two records record1 and record2 in a list.
    @param input_protocol: A number from Constants indicating the input protocol.
    @param layout: A number from Constants indicating the type of layout.
    @param methylation_counter: A dictionary for counting methylation context information.
    @param genome_dictionary: A dictionary of the reference genome.  Set to None for multiprocessing.
    @param genome_pipe: A pipe to ask for genome reference information for multiprocessing.
    @return A tuple object for writing a record, the number of multireads uniquified, and the alignment direction information.
    """
    ambiguous_uniquified = 0
    direction_id = None
    if empty: #An unmapped SAM record will be written if the FASTQ record could not be found in the SAM file.  This shouldn't happen, but it probably will.
        write_record = (Constants.SAM_WRITE_UNMAPPED, [[BisPin_util.changeFastqRecordToSAMRecord(record1, seq_id = alignement_key)]])
        if record2 != None:
            write_pair = BisPin_util.changeFastqRecordToSAMRecord(record2, seq_id = alignement_key)
            write_record[1][0].append(write_pair)
    else:
        #Choose the SAM records with the highest alignment scores.
        for key in meilleur_alignement_dictionnaire:
                meilleur_alignement_dictionnaire[key] = choisirAlignement(meilleur_alignement_dictionnaire[key])
        meilleur_alignements = choisirEntreDictionnaires(meilleur_alignement_dictionnaire, 
                                                         record1, record2, alignement_key, 
                                                         cutoff_bs, cutoff_r, input_protocol)
        if meilleur_alignements[0] == Constants.SAM_WRITE_ALIGNED: 
            if rescore_matrix[0] != None:
                len_ma1 = len(meilleur_alignements[1])
            #Process the best SAM records adding the methylation information, updating the FLAG, etc.
            meilleur_alignements = traiterAlignement(meilleur_alignements, record, 
                                                     input_protocol, layout, methylation_counter, 
                                                     rescore_matrix, genome_dictionary, genome_pipe = genome_pipe)
            len_ma2 = len(meilleur_alignements)
            if rescore_matrix[0] != None and len_ma2 == 1 and len_ma1 > 1: #Detects if a multiread was uniquified.
                ambiguous_uniquified += 1
            if len_ma2 > 1:  #What about hairpin?
                write_record = (Constants.SAM_WRITE_AMBIGUOUS, meilleur_alignements)
            else:
                read_direction = ""
                for record in meilleur_alignements[0]:
                    read_direction += "_" + record[Constants.SAM_KEY_READ]
                genome_direction = meilleur_alignements[0][0][Constants.SAM_KEY_GENOME]
                direction_id = read_direction[1:] + "_" + genome_direction
                write_record = (Constants.SAM_WRITE_UNIQUE, meilleur_alignements)
        elif meilleur_alignements[0] == Constants.SAM_WRITE_UNMAPPED:
            write_record = (Constants.SAM_WRITE_UNMAPPED, meilleur_alignements[1])
        elif meilleur_alignements[0] == Constants.SAM_WRITE_FILTERED:
            for paired_record in meilleur_alignements[1]:
                for single_record in paired_record:
                    for key in list(single_record.keys()):
                        if not key in Constants.SAM_KEYS_TO_KEEP:
                            del single_record[key]
            write_record = (Constants.SAM_WRITE_FILTERED, meilleur_alignements[1])
    return (write_record, ambiguous_uniquified, direction_id)


def makeSamFileTypeDictionary(sam_file_list):
    """
    Create a dictionary of SAM file iterator objects.  The SAM file must follow a specific naming convention.
    @param sam_file_list: A list of strings representing the location of raw BFAST SAM files.
    @return: A dictionary of the SAM file iterator objects and a dictionary giving the conversion information for the SAM files. 
    """
    filetypes = Constants.CONV_TYPES
    sam_iterator_dict = {}
    sam_conversion_dict = {}
    onlyfiles = sam_file_list
    for a_file in onlyfiles:
        split_a_file = a_file.split(".")
        for i in range(len(split_a_file)-1, -1, -1):
            filetypes_equal = list(filter((lambda x: x == split_a_file[i]), filetypes))
            if filetypes_equal != [] and filetypes_equal != '':
                s_type = filetypes_equal[0]
                sam_iterator_dict[s_type] = SeqIterator.SeqIterator(a_file, file_type=Constants.SAM)
                sam_conversion_dict[s_type] = s_type.split("_")
                break
    return (sam_iterator_dict, sam_conversion_dict)


def choisirAlignement(brut_alignements):
    """
    Chooses the highest scoring alignment(s) from a list of SAM records.  Returns a list of lists of grouped SAM records.
    This is the format used throughout the post processing module to represent SAM records after this function has been called. 
    @param brut_alignements: a list of SAM records from the same SAM file. (Each SAM record is a dictionary)
    @return: a tuple of the best value and an iterable of lists of SAM records that represent the best alignment(s)
    Each list of SAM records in the iterable can either be a 2-list for paired end data or a 1-list for single end data.
    e.g. [[pair1, pair2], [single1]]
    If the alignments are all unmapped, then return None
    """
    unmapped_bits = list(map((lambda record: isUnmapped(record[Constants.SAM_KEY_FLAG])), 
                         brut_alignements))
    all_unmapped = reduce((lambda x, y: x and y), unmapped_bits)
    if all_unmapped: #If all the records are unmapped, return None
        return None
    else:
        grouped_records = []
        firstSegments = []
        lastSegments = []
        #Search through each record and group them according to whether they are paired or not.
        for record in brut_alignements:
            flag = record[Constants.SAM_KEY_FLAG]
            if hasMultipleSegments(flag) and isFirstSegment(flag):
                firstSegments.append(record)
            elif hasMultipleSegments(flag) and isLastSegment(flag):
                lastSegments.append(record)
            else:
                grouped_records.append([record])
        grouped_records += [list(group) for group in zip(firstSegments, lastSegments)]
        valeurs_et_records = []
        #Get the average alignment score for each group.
        for group in grouped_records:
            group_sum  = reduce(lambda x, y: x + float(y[Constants.SAM_KEY_ALIGNMENT_SCORE]), group, 0.0)
            group_avg = group_sum / (len(group) + 0.0)
            valeurs_et_records.append((group_avg, group))
    #Find the maximum alignment score and return only those grouped records.
    valeur_max = max(valeurs_et_records, key = (lambda x : x[0]))[0]
    valeurs_et_records_max = list(filter((lambda x: x[0] == valeur_max), valeurs_et_records))
    retvalue = (valeur_max, list(map((lambda x: x[1]), valeurs_et_records_max)))
    return retvalue
        
        
def isUnmapped(flag):
    """Checks if the SAM flag indicates an unmapped read."""
    return ((int(flag) >> 2) % 2) == 1

def isReversed(flag):
    """Checks if the SAM flag indicates a reversed read."""
    return ((int(flag) >> 4) % 2) == 1

def setFilterBit(flag):
    """Sets the filter bit for a SAM flag."""
    return int(flag) | 512
    
def hasMultipleSegments(flag):
    """Checks if the SAM flag indicates multiple segments."""
    return (int(flag) % 2) == 1

def isFirstSegment(flag):
    """Checks if the SAM flag indicates the first segment."""
    return ((int(flag) >> 6) % 2) == 1

def isLastSegment(flag):
    """Checks if the SAM flag indicates the last segment."""
    return ((int(flag) >> 7) % 2) == 1


def choisirEntreDictionnaires(meilleur_alignement_dictionnaire, fastq_record1, fastq_record2, 
                              alignement_key, cutoff_bs, cutoff_r, input_protocol):
    """
    Chooses the best scoring SAM records from different SAM files representing different conversion protocols
    @param: meilleur_alignement_dictionnaire: A dictionary where the key represents a SAM file and the value is the best SAM records
    for that key corresponding to one FASTQ record.
    @param: fastq_record1: The first FASTQ record.
    @param: fastq_record2: None for single end data.  The second FASTQ record for paired end data.
    @param: alignement_key: The sequence id for the FASTQ record.  This should be identical to the SAM records.
    @param: cutoff_bs: A float for filtering low quality reads.
    @param: cutoff_r: A float for filtering low quality hairpin recovered reads.
    @param: input_protocol: A number from Constants indicating the input protocol.
    @return: a tuple where the first position indicates the type of read and the second
    position is an iterable of SAM records.  (Each SAM record is a dictionary) 
    """
    alignement_liste = [meilleur_alignement_dictionnaire[key] 
                            for key in meilleur_alignement_dictionnaire 
                            if meilleur_alignement_dictionnaire[key] != None]
    if alignement_liste == []: #If the dictionary is empty, then the reads are unmapped.
        grouped_record = [BisPin_util.changeFastqRecordToSAMRecord(fastq_record1, seq_id = alignement_key)]
        if fastq_record2 != None:
            grouped_record.append(BisPin_util.changeFastqRecordToSAMRecord(fastq_record2, seq_id = alignement_key))
        write_record =  (Constants.SAM_WRITE_UNMAPPED, [grouped_record])
        return write_record
    #Find the maximum score for the records.
    len_seq = float(len(fastq_record1[1]))
    if fastq_record2 != None:
        len_seq += len(fastq_record2[1])
        len_seq /= 2.0
    valeur_max = max(alignement_liste, key=(lambda x: x[0]))[0]
    meilleur_alignement_liste = []
    records_added = 0
    conversion_type = None
    #Choose records that have the highest alignment score.
    for m in alignement_liste:
        if m[0] == valeur_max:
            if conversion_type == None or conversion_type == Constants.CONV_R:
                conversion_type = m[1][0][0][Constants.SAM_KEY_READ]
            meilleur_alignement_liste += m[1]
            records_added += 1
    #Report these SAM records as filtered if the maximum alignment score is too low.  Otherwise report aligned.
    if (conversion_type == Constants.CONV_R
         and float(valeur_max) / (len_seq + 0.0) <= cutoff_r) or (
                                                                  conversion_type != Constants.CONV_R and 
                                                                  float(valeur_max) / (len_seq + 0.0) <= cutoff_bs):
        return (Constants.SAM_WRITE_FILTERED, changerFiltre(meilleur_alignement_liste))
    else:        
        return (Constants.SAM_WRITE_ALIGNED, meilleur_alignement_liste)
       
            
def changerFiltre(sam_record_liste):
    """
    Transforms the SAM records in the input list into filtered (low quality alignment) SAM records
    This is done by setting a bit in the FLAG field. 
    @param sam_record_liste: a list of of lists of grouped SAM records.  (Each SAM record is a dictionary)
    @return: a list of of lists of grouped SAM records that represent a filtered set of records
    """
    for record in sam_record_liste:
            for subrecord in record:
                subrecord[Constants.SAM_KEY_FLAG] = str(int(subrecord[Constants.SAM_KEY_FLAG]) 
                                                        | int(Constants.SAM_VALUE_FILTERED))
    return sam_record_liste


def traiterAlignement(meilleur_alignement, fastq_record_group, input_protocol, layout, methylation_count, 
                      rescore_matrix, genome_dictionary, genome_pipe = None):
    """
    Processes alignments corresponding to one FASTQ record by modifying the MD tag, updating the flag, 
    computing the methylation calling string, and other things.
    @param meilleur_alignement: SAM alignement records to be processed that correspond to the same FASTQ record
    @param fastq_record_group: A list of FASTQ records.  There are two records for paired end data.
    @param input_protocol: A number from Constants indicating the input protocol.
    @param layout: A number from Constants indicating the type of layout.
    @param methylation_count: A dictionary for counting methylation context information.
    @param rescore_matrix: A 2-tuple of dictionaries representing the rescoring function for multireads.  If both elements are set to None, no rescoring will be done.
    @param genome_dictionary: A dictionary of the reference genome.  Set to None for multiprocessing.
    @param genome_pipe: A pipe to ask for genome reference information for multiprocessing.
    @return: the processed SAM records as a list of lists of grouped SAM dictionary records
    """
    len_meilleur_alignement = len(meilleur_alignement[1])
    is_ambiguous = (len_meilleur_alignement > 1)
    #Look at each grouped SAM record
    for paired_record in meilleur_alignement[1]:
        #Look at the SAM record and its corresponding FASTQ record from the input
        for subrecord, fastq_record in zip(paired_record, fastq_record_group):
            #Delete unnecessary fields.
            for key in list(subrecord.keys()):
                if not key in Constants.SAM_KEYS_TO_KEEP:
                    del subrecord[key]
            #Extract the CIGAR string and the MD tag for computing and modifying the alignment strings
            token_cigar = BisPin_util.tokenizeCigar(subrecord[Constants.SAM_KEY_CIGAR])
            token_md = BisPin_util.tokenizeMDtag(subrecord[Constants.SAM_KEY_MD])
            read_conversion = subrecord[Constants.SAM_KEY_READ]
            genome_conversion = subrecord[Constants.SAM_KEY_GENOME]
            #There is a bug in BFAST's SAM flag that reports the incorrect orientation
            #for the second segment for paired end reads.
            #These flags are for the GA converted genome
            #115 should be 83 RC
            #179 should be 147 not RC
            #These flags are for the CT converted genome
            #65 should be 99 not RC
            #129 should be 147 RC
#             if isReversed(subrecord[Constants.SAM_KEY_FLAG]):
#                 direction = Constants.CONV_REVERSE
            if read_conversion != genome_conversion:
                direction = Constants.CONV_REVERSE
            else:
                direction = Constants.CONV_FORWARD
            fastq_seq = fastq_record[1]
            len_fastq_record = len(fastq_seq)
            try: #Get the reference sequence and little extra (cdr) to determine context.
                reference_sequence, ref_cdr = getReferenceSequence(subrecord[Constants.SAM_KEY_RNAME].split()[0], 
                                                      int(subrecord[Constants.SAM_KEY_POS]), 
                                                      len_fastq_record, token_cigar, 
                                                      read_conversion, genome_conversion, direction,
                                                      genome_dictionary, genome_pipe = genome_pipe)
            except TypeError:
                sys.stderr.write(logstr + "A problem occured getting the reference sequence\n")
                sys.stderr.write(logstr + "fastq_record\t%s\n" % str(fastq_record))
                raise
            if direction == Constants.CONV_REVERSE:
                fastq_seq = BisPin_util.reverseComplement(fastq_seq)
            #Compute the mismatch, indel, and deletion list (MI_string) and the new MD tag and the hamming distance. 
            new_md, MI_string, hamming_distance = makeMismatchString(fastq_seq, 
                                        subrecord, 
                                        reference_sequence,
                                        token_cigar, 
                                        token_md)
            #Rescore the alignment with the rescore_matrix
            if is_ambiguous and rescore_matrix[0] != None:
                subrecord[Constants.SAM_KEY_RESCORE] = str(rescoreAlignment(
                                                                                    fastq_seq, MI_string, direction, 
                                                                                    rescore_matrix, read_conversion, genome_conversion))
            #Create the methylation calling string.
            if read_conversion == Constants.CONV_R:
                subrecord[Constants.SAM_KEY_RECOVERED] = subrecord[Constants.SAM_KEY_SEQ]
                (methylation_CT, methylation_GA) = hairpinRecoveredMethylationCall(fastq_record_group[0][1],
                                                                                   fastq_record_group[1][1],
                                                                                   subrecord[Constants.SAM_KEY_RECOVERED])
                subrecord[Constants.SAM_KEY_METHYLATION] = methylation_CT
                subrecord[Constants.SAM_KEY_METHYLATION_GA] = methylation_GA
            else:
                methyl_string = makeMethylationCallString(fastq_seq, 
                                                      ref_cdr, 
                                                      MI_string,
                                                      direction,
                                                      read_conversion)
                subrecord[Constants.SAM_KEY_METHYLATION] = methyl_string
            subrecord[Constants.SAM_KEY_MD] = new_md
            subrecord[Constants.SAM_KEY_SEQ] = fastq_seq
            subrecord[Constants.SAM_KEY_PROGRAM] = Constants.SAM_VALUE_PROGRAM
            subrecord[Constants.SAM_KEY_DISTANCE] = str(hamming_distance)
    if is_ambiguous and rescore_matrix[0] != None: # selectionnez le meilleur alignement pour multireads
        valeurs_et_records = list(map(calculateScore, meilleur_alignement[1]))
        valeur_max = max(valeurs_et_records, key = (lambda x : x[0]))[0]
        valeurs_et_records_max = list(filter((lambda x: x[0] == valeur_max), valeurs_et_records))
        meilleur_alignement_choice = list(map((lambda x: x[1]), valeurs_et_records_max))
    else:
        meilleur_alignement_choice = meilleur_alignement[1]
    if is_ambiguous or rescore_matrix[0] != None: #Adjust the bit on the flag to indicate the presence of secondary alignments
        primary_locations = [i for i in range(len(meilleur_alignement_choice)) if reduce(lambda x, y: x and y ,  
                                                                                         [((Constants.SAM_VALUE_SECONDARY &
                                                                                                int(single_record[Constants.SAM_KEY_FLAG])) >> 8) ==
                                                                                              0 for single_record in meilleur_alignement_choice[i]]) ]
        num_primaries = len(primary_locations)
        if num_primaries == 0: # No primary sequences.  Set a sequence to be the primary one.
            for subrecord in meilleur_alignement_choice[0]:
                subrecord[Constants.SAM_KEY_FLAG] = str(int(subrecord[Constants.SAM_KEY_FLAG]) ^ Constants.SAM_VALUE_SECONDARY)
        elif num_primaries > 1: # Change other sequences to be secondary sequences.
            for i in primary_locations[1:num_primaries]:
                for subrecord in meilleur_alignement_choice[i]:
                    subrecord[Constants.SAM_KEY_FLAG] = str(int(subrecord[Constants.SAM_KEY_FLAG]) | Constants.SAM_VALUE_SECONDARY)
    if meilleur_alignement[0] == Constants.SAM_WRITE_ALIGNED: #Count the methylation marks
        for alignement_record in meilleur_alignement_choice:
            for subrecord in alignement_record:
                    for char in subrecord[Constants.SAM_KEY_METHYLATION]:
                        methylation_count[char] = methylation_count.get(char, 0) + 1
    return meilleur_alignement_choice


def calculateScore(record, key = Constants.SAM_KEY_RESCORE):
    """
    Return the average alignment score for grouped SAM dictionary records.
    @param record: A grouping (list) of SAM dictionary records
    @param  key: Governs which alignment score to use.  Defaults to the rescored score.
    @return: The average alignment score for grouped SAM records.
    """
    if isinstance(record, list):
        score = sum([float(x[key]) for x in record])
        score /= (len(record) + 0.0)
        return (score, record)
    else:
        return (float(record[key]), record)


def rescoreAlignment(fastq_seq, MI_string, direction, score_matrix, read_conversion, genome_conversion):
    """"
    Rescores an alignment for uniquifying multireads
    @param fastq_seq: The FASTQ sequence
    @param MI_string: The MI_string which represents the alignment.  Includes mismatches and indels.
    @param direction: Whether the alignment is forward or reverse.
    @param score_matrix: A 2-tuple of dictionaries representing the rescoring matrix.
    @param read_conversion: (Not used.)
    @param genome_conversion: A string representing the genome conversion for the alignment.
    @return Returns a new score for an alignment based on the my_score_matrix
    """
    if score_matrix[1] == None or genome_conversion == Constants.CONV_CT:
        my_score_matrix = score_matrix[0]
    else:
        my_score_matrix = score_matrix[1] 
    score = sum(  #Add the contributions of deletions.
                map((lambda x: my_score_matrix[Constants.OPEN_DEL] + my_score_matrix[Constants.EXT_DEL]*(len(x) - 1)) , 
                    list(filter((lambda x: x.startswith("^")), MI_string))))
    #Filter out the deletions so that insertions and mismatches can be counted
    MI_string = list(filter((lambda x: not x.startswith("^")), MI_string))
    openInsert = False
    for i in range(len(MI_string)):
        MI_char = MI_string[i].upper()
        if MI_char == Constants.CIGAR_I: #Calculate the gap penalty.
            score += my_score_matrix[Constants.EXT_INS]
            if not openInsert:
                openInsert = True
                score += my_score_matrix[Constants.OPEN_INS]
        else: #Calculate the match or mismatch score
            openInsert = False
            if MI_char == Constants.CIGAR_M: #Match
                score += my_score_matrix[fastq_seq[i].upper()][fastq_seq[i].upper()]
            elif genome_conversion == Constants.CONV_GA and MI_char == 'G' and fastq_seq[i].upper() == 'A': #GA methylation
                score += my_score_matrix["GA"]
            elif genome_conversion == Constants.CONV_CT and MI_char == 'C' and fastq_seq[i].upper() == 'T': #CT methylation
                score += my_score_matrix["CT"]
            else: #Mismatch
                try:
                    score += my_score_matrix[MI_char][fastq_seq[i].upper()]
                except KeyError:
                    sys.stderr.write("%sThe rescoring function did not recognize a character at %d.  Using 'N' instead.  MI_string: %s  fastq_seq: %s\n" % (logstr, i, str(MI_string), fastq_seq))
                    score += my_score_matrix["N"][fastq_seq[i].upper()]
    return score
            

def getAlignmentDirection(read_direction, genome_direction):
    """
    Returns whether the alignment should be a forward or a reverse alignment.
    This function uses string constants found in the Constants module.
    @param read_direction: The conversion state of the reads.
    @param genome_direction: The conversion state of the reference genome.
    @return: A string representing the alignment direction 
    """
    if (read_direction == Constants.CONV_CT and genome_direction == Constants.CONV_CT) or \
        (read_direction == Constants.CONV_GA and genome_direction == Constants.CONV_GA) or \
        (read_direction == Constants.CONV_R):
        return Constants.CONV_FORWARD
    else:
        return Constants.CONV_REVERSE
    
    
# 
# def getReferenceSequence(key, position, seq_len, token_cigar, direction, genome_dictionary, genome_pipe = None):
#     """
#     Gets the reference sequence for a SAM record.  This is so that the MD tag can be updated.
#     
#     @param key: A key for the sequence (the id for the fasta file record)
#     @param position: The starting position of the string starting at the string at the key location
#     in the genome dictionary
#     @param seq_len: The length of the sequence in the SAM record
#     @param token_cigar: The tokenized cigar string from the unprocessed SAM record
#     @param direction: A string determining whether the alignment corresponding to the requested 
#     reference sequence is forward or reverse. 
#     @param genome_dictionary: A dictionary representing the reference fasta file
#     @param genome_pipe: A pipe to send and receive genome pipe information if access to the genome
#     is controlled by another process. (Not used.)
#     @return: The reference sequence that corresponds to the SAM record.
#     """
#     position = int(position)
#     len_insertions = 0
#     len_deletions = 0
#     for token in token_cigar:
#         if token[1] == Constants.CIGAR_D:
#             len_deletions += token[0]
#         elif token[1] == Constants.CIGAR_I:
#             len_insertions += token[0]
#     reference_sequence_length = seq_len - len_insertions + len_deletions
#     position -= 1
#     if direction == Constants.CONV_FORWARD:
#         index1 = position
#         index2 = position + reference_sequence_length + 2
#     elif  direction == Constants.CONV_REVERSE and position >= 2:
#         index1 = position - 2
#         index2 = position + reference_sequence_length
#     else:
#         index1 = position
#         index2 = position + reference_sequence_length
#     if genome_pipe == None:
#         ref_string =  genome_dictionary[key][index1: index2]
#     else:
#         genome_pipe.send((key, index1, index2))
#         ref_string = genome_pipe.recv()
#     if direction == Constants.CONV_FORWARD:
#         ref_string_car = ref_string[0:reference_sequence_length]
#         ref_string_cdr = ref_string[reference_sequence_length: reference_sequence_length + 2]
#         if len(ref_string_cdr) < 2:
#             ref_string_cdr = 'XX'
#     elif direction == Constants.CONV_REVERSE and position >=2:
#         ref_string_car = ref_string[2:len(ref_string)]
#         ref_string_cdr = "".join(reversed(ref_string[0:2]))
#     else:
#         ref_string_car = copy.copy(ref_string)
#         ref_string_cdr = 'XX'
#     return (ref_string_car.upper(), ref_string_cdr.upper())


def getReferenceSequence(key, position, seq_len, token_cigar, read_conversion, genome_conversion, 
                         direction, genome_dictionary, genome_pipe = None):
    """
    Gets the reference sequence for a SAM record.  This is so that the MD tag can be updated.
    
    @param key: A key for the sequence (the id for the fasta file record)
    @param position: The starting position of the string starting at the string at the key location
    in the genome dictionary
    @param seq_len: The length of the sequence in the SAM record
    @param token_cigar: The tokenized cigar string from the unprocessed SAM record
    @param direction: A string determining whether the alignment corresponding to the requested 
    reference sequence is forward or reverse. 
    @param genome_dictionary: A dictionary representing the reference fasta file
    @param genome_pipe: A pipe to send and receive genome pipe information if access to the genome
    is controlled by another process. (Not used.)
    @return: The reference sequence that corresponds to the SAM record.
    """
    position = int(position)
    len_insertions = 0
    len_deletions = 0
    for token in token_cigar:
        if token[1] == Constants.CIGAR_D:
            len_deletions += token[0]
        elif token[1] == Constants.CIGAR_I:
            len_insertions += token[0]
    reference_sequence_length = seq_len - len_insertions + len_deletions
    position -= 1
    if genome_pipe == None:
        ref_string_car = genome_dictionary[key][position: position + reference_sequence_length]
    if (read_conversion == Constants.CONV_CT and genome_conversion == Constants.CONV_CT) or \
        (read_conversion == Constants.CONV_GA and genome_conversion == Constants.CONV_CT) or \
        (read_conversion == Constants.CONV_R and direction == Constants.CONV_FORWARD): #Get the upstream cdr
        if genome_pipe != None:
            genome_pipe.send((key, position, position + reference_sequence_length + 2))
            ref_string = genome_pipe.recv()
            ref_string_car = ref_string[0:reference_sequence_length]
            ref_string_cdr = ref_string[reference_sequence_length:reference_sequence_length + 2]
        else:
            ref_string_cdr = genome_dictionary[key][position + reference_sequence_length : position + reference_sequence_length + 2]
        if len(ref_string_cdr) < 2: #The upstream cdr is near the end of the reference sequence.
            ref_string_cdr = 'NN'
    elif position >= 2 and ((read_conversion == Constants.CONV_CT and genome_conversion == Constants.CONV_GA) or \
        (read_conversion == Constants.CONV_GA and genome_conversion == Constants.CONV_GA) or \
        (read_conversion == Constants.CONV_R and direction == Constants.CONV_REVERSE)): #Get the downstream cdr
        if genome_pipe != None:
            genome_pipe.send((key, position - 2, position + reference_sequence_length))
            ref_string = genome_pipe.recv()
            ref_string_car = ref_string[2 : 2 + reference_sequence_length]
            ref_string_cdr = ref_string[0 : 2]
        else:
            ref_string_cdr = genome_dictionary[key][position - 2 : position]
        ref_string_cdr = ref_string_cdr[::-1]
    else: # The downstream cdr is at the beginning of the reference sequence.
        if genome_pipe != None:
            genome_pipe.send((key, position, position + reference_sequence_length))
            ref_string_car = genome_pipe.recv()
        ref_string_cdr = 'NN'
    return (ref_string_car.upper(), ref_string_cdr.upper())


def makeReferenceGenomeDictionary(reference_genome_file_location):
    """
    Makes a dictionary for the reference genome indexed by contig id.
    @param reference_genome_file_location: A string representing the location of the reference genome file.
    @return: a dictionary object indexed by contig id where the value is the reference sequence, 
    A list of ids (keys) and the length of each contig, and the time it took to load the genome.
    """
    now = datetime.datetime.now()
    genome_iterator = SeqIterator.SeqIterator(reference_genome_file_location, file_type=Constants.FASTA)
    id_length_list = []
    genome_dictionary = {}
    for record in genome_iterator:
        key = record[0].split()[0]
        genome_dictionary[key] = record[1]
        id_length_list.append((key, str(len(record[1]))))
    later = datetime.datetime.now()
    return (genome_dictionary, id_length_list, later - now)


def checkDeletionInsertionProblem(token_cigar, token_md):
    """
    This functions checks that the length of any deletion in the token md list is the 
    same as the length of the deletion in the cigar list.  The MD list is corrected if 
    there is a discrepancy.  This hacks a solution to a bug encountered in BFAST-Gap.
    @param token_cigar: The token list representing the cigar string 
    @param param: The The token list representing the cigar string
    @return A token cigar list followed by a corrected token md list
    """
    differences = []
    start = 0
    for i in range(len(token_cigar)):
        if token_cigar[i][1] == 'D':
            for j in range(start, len(token_md)):
                if str(token_md[j]).startswith('^') and int(token_cigar[i][0]) != len(token_md[j]) - 1:
                    differences.append((int(token_cigar[i][0]), j))
                elif str(token_md[j]).startswith('^'):
                    start = j + 1
                    break
    for diff in differences:
        i,j = diff
        token_md[j] = token_md[j][0:i + 1]
        if j+1 >= len(token_md):
            token_md.append(1)
        elif isinstance(token_md[j+1], int):
            token_md[j+1] += 1
        else:
            token_md.insert(j+1, 1)
    return (token_cigar, token_md)


def makeMismatchString(bisulfite_sequence, alignement_record, reference_sequence, cigar, md, printOut = False):
    """
    @param bisulfite_sequence: is the unconverted short read
    @param alignement_record: a SAM record dictionary for making the mismatch (alignment) string.
    @param reference_sequence: the unconverted reference sequence that matches the location for the alignment
    @param cigar: is the tokenized cigar string from the sam output
    @param md: is the tokenized MD tag string from the sam output
    @param printOut: A boolean that controls whether debug code is printed out.
    @return: the md_tag string for the alignment between the bisulfite_sequence and the converted_sequence
    followed by the MI_string (list) that is used to create the md_tag, and the hamming distance of the alignment.
    """
    converted_sequence = alignement_record[Constants.SAM_KEY_SEQ]
    if len(converted_sequence) != len(bisulfite_sequence):
        sys.stderr.write("%smakeMismatchString: Lengths of converted sequence and bisulfite sequence do not match.\nAlignment Record:\t%s\nInput FASTQ Record:\t%s\n" %
        (logstr, str(alignement_record), str(bisulfite_sequence)))
    token_cigar = BisPin_util.tokenizeCigar(cigar) if isinstance(cigar, str) else cigar
    token_md = BisPin_util.tokenizeMDtag(md) if isinstance(md, str) else md
    token_cigar, token_md = checkDeletionInsertionProblem(token_cigar, token_md)
    deletions_list = BisPin_util.findDeletions(token_cigar, token_md, reference_sequence)
    if printOut:
        sys.stderr.write(logstr + deletions_list + "\n")
    MI_string = []
    #Create the initial MI string consisting of Ms and Is
    for token in token_cigar: 
        if token[1] == 'M' or token[1] == 'I':
            MI_string.extend(token[0] * token[1])
    if printOut:
        sys.stderr.write(logstr + MI_string + "\n")
    #Add deletions to the MI_string list
    if len(deletions_list) != 0: 
        new_MI_string = []
        counter_deletes = 0
        counter_md = 0
        for i in range(len(MI_string)):
            if (MI_string[i] != 'I'):
                counter_md += 1
            if len(deletions_list) > counter_deletes and deletions_list[counter_deletes][0] == counter_md:
                new_MI_string.append(deletions_list[counter_deletes][1])
                counter_deletes += 1
            new_MI_string.append(MI_string[i])
        MI_string = new_MI_string
        if printOut:
            sys.stderr.write(logstr + MI_string + "\n")
    bisulfite_index = 0
    reference_index = 0
    hamming_distance = 0
    #Check for mismatches between the reference and the bisulfite sequence
    for index in range(len(MI_string)): 
        if MI_string[index] == Constants.CIGAR_I:
            bisulfite_index += 1
            hamming_distance += 1
        elif MI_string[index].startswith('^'):
            del_length = len(MI_string[index]) - 1
            reference_index += del_length
            hamming_distance += del_length
        else:
            try:
                bisulfite_char = bisulfite_sequence[bisulfite_index]
            except IndexError:
                error_string = "****BSSEQ_ERROR*****\nbisulfite_sequence:\t%s\nreference_sequence:\t%s\nindex:\t%s\ncigar:\t%s\nmd:\t%s\nMI_string:\t%s\nlen(MI_string):\t%s\nalignement_record:\t%s\nbisulfite_index:\t%s\n****ERROR*****\n" % (str(bisulfite_sequence), str(reference_sequence), str(index), str(cigar), str(md), str(MI_string), str(len(MI_string)), str(alignement_record), str(bisulfite_index))
                sys.stderr.write(logstr + error_string + "\n")
                sys.stderr.flush()
                raise
            try:
                reference_char = reference_sequence[reference_index]
            except IndexError:
                error_string = "****REFSEQ_ERROR*****\nbisulfite_sequence:\t%s\nreference_sequence:\t%s\nindex:\t%s\ncigar:\t%s\nmd:\t%s\nMI_string:\t%s\nlen(MI_string):\t%s\nalignement_record:\t%s\nreference_index:\t%s\n****ERROR*****\n" % (str(bisulfite_sequence), str(reference_sequence), str(index), str(cigar), str(md), str(MI_string), str(len(MI_string)), str(alignement_record), str(reference_index))
                sys.stderr.write(logstr + error_string + "\n")
                sys.stderr.flush()
                raise
            if bisulfite_char != reference_char:
                MI_string[index] = reference_char
                hamming_distance += 1
            else:
                MI_string[index] = Constants.CIGAR_M
            bisulfite_index += 1
            reference_index += 1
    if printOut:
        sys.stderr.write(logstr + MI_string + "\n")
    new_md_string = ""
    num_mi = 0
    last_char = ""
    #Generate the new MD:Z tag string
    for i in range(len(MI_string)): 
        MI_char = MI_string[i]
        if last_char == "" and MI_char != 'M' and MI_char != 'I':
            new_md_string += MI_char
        elif (MI_char == 'M'): 
            num_mi += 1
        elif MI_char != 'M' and MI_char != 'I' and last_char == 'M':
            new_md_string += str(num_mi)
            num_mi = 0
            new_md_string += MI_char
        elif last_char != 'M' and last_char != 'I' and MI_char != 'M' and MI_char != 'I':
            new_md_string += str(0)
            new_md_string += MI_char
        elif last_char == 'I' and MI_char != 'M' and MI_char != 'I':
            if (num_mi > 0):
                new_md_string += str(num_mi)
                num_mi = 0
            new_md_string += MI_char
        last_char = MI_char
    if num_mi > 0:
        new_md_string += str(num_mi)
        num_mi = 0
    if len(new_md_string) > 0:
        if not new_md_string[0].isdigit():
            new_md_string = '0' + new_md_string
        if not new_md_string[-1].isdigit():
            new_md_string = new_md_string + '0'
    return new_md_string, MI_string, hamming_distance


def hairpinRecoveredMethylationCall(bisulfite_sequence_CT, bisulfite_sequence_GA, recovered_original):
    """
    Gives methylcation calling string for recovered hairpin sequences.  Gives both the 
    CT converted and GA converted methylation calling strings.
    This uses the recovered original sequence rather than the reference genome so that the 
    methylation call is taken from the organism that generated the sequence.
    
    @see: makeMethylationCallString documentation for the definition of methylation contexts 
    @param bisulfite_sequence_CT: The C to T converted bisulfite sequence for hairpin data
    @param bisulfite_sequence_GA: The G to A converted bisulfite sequence for hairpin data
    @param recovered_original: The recovered untreated original strand 
    @return The methylation calling string corresponding to the CT converted strand and 
    the methylation calling string corresponding to the GA converted strand
    """
    len_diff = abs(len(bisulfite_sequence_CT) - len(bisulfite_sequence_GA))
    extra = "N" * (len_diff + 2) #Add extra 'X' characters so that characters aren't run out
    min_len = min(len(bisulfite_sequence_CT), len(bisulfite_sequence_GA))
    methylation_list = []
    recovered_original_X = recovered_original + extra
    for i in range(min_len): # Find the methylation calling string for the C to T strand. 
        first_char = recovered_original_X[i + 1]
        second_char = recovered_original_X[i + 2]
        if bisulfite_sequence_CT[i] == 'C' and bisulfite_sequence_GA[i] == 'C':
            if first_char == 'G':
                methylation_list.append('Z')
            elif second_char == 'G':
                methylation_list.append('X')
            elif first_char != 'N' and second_char != 'N':
                methylation_list.append('H')
            else:
                methylation_list.append('U')
        elif bisulfite_sequence_CT[i] == 'T' and bisulfite_sequence_GA[i] == 'C':
            if first_char == 'G':
                methylation_list.append('z')
            elif second_char == 'G':
                methylation_list.append('x')
            elif first_char != 'N' and second_char != 'N':
                methylation_list.append('h')
            else:
                methylation_list.append('u')
        else:
            methylation_list.append('.')
    methylation_CT = "".join(methylation_list)
    #The G to A converted strand needs to be reversed so that the methylation calling context goes form 5' to 3' 
    bisulfite_sequence_GA_reverse = bisulfite_sequence_GA[0:min_len][::-1]
    recovered_original_reversed = recovered_original[0:min_len][::-1] + extra
    bisulfite_sequence_CT_reverse = bisulfite_sequence_CT[0:min_len][::-1]
    methylation_list = []
    for i in range(min_len): #Compute the methylation calling string for the G to A strand.
        first_char = recovered_original_reversed[i + 1]
        second_char = recovered_original_reversed[i + 2]
        if bisulfite_sequence_CT_reverse[i] == 'G' and bisulfite_sequence_GA_reverse[i] == 'G':
            if first_char == 'G':
                methylation_list.append('Z')
            elif second_char == 'G':
                methylation_list.append('X')
            elif first_char != 'N' and second_char != 'N':
                methylation_list.append('H')
            else:
                methylation_list.append('U')
        elif bisulfite_sequence_CT_reverse[i] == 'G' and bisulfite_sequence_GA_reverse[i] == 'A':
            if first_char == 'G':
                methylation_list.append('z')
            elif second_char == 'G':
                methylation_list.append('x')
            elif first_char != 'N' and second_char != 'N':
                methylation_list.append('h')
            else:
                methylation_list.append('u')
        else:
            methylation_list.append('.')
    methylation_GA = "".join(methylation_list[::-1])
    return (methylation_CT, methylation_GA)
    
    
def makeMethylationCallString(bisulfite_sequence, ref_cdr, MI_string, direction, read_conversion):
    """
    Makes a methylation calling string.  The following text gives the format used.
    Z - methylated C in the CpG context
    z - unmethylated C in the CpG context
    X - methylated C in the CHG context
    x - unmethylated C in the CHG context
    H - methylated C in the CHH context
    h - unmethylated C in the CHH context
    U - methylated C in unknown CN or CHN context 
    u - unmethylated C in unknown CN or CHN context
    . - not a cytosine 
    Bismark uses the reference sequence to determine the context, so this code also uses the reference sequence.
    Using the bisulfite sequence called 70 percent CpG methylated on a simulation for 80 percent CpG methylated.
    Using the bisulfite sequence might make sense as the context is called relative to the organism and not the reference;
    however, what is to be done about the bisulfite conversion, etc.?
    @param bisulfite_sequence: An unconverted bisulfite sequence
    @param ref_cdr: The endpart of the reference sequence to be concatenated onto the bisulfite_sequence 
    so that it is not over run.
    @param MI_string: A list representing an alignment.  Includes info on mismatchtes and indels.
    This can be used to get the base from the reference genome.
    @param direction: The direction of the alignment (forward or reverse).  Uses constants from Constants.
    @return: The methylation call string.
    """
    MI_string = list(filter((lambda x: not x.startswith("^")), MI_string))
    methylation_list = []
    len_MI_string = len(MI_string)
    #Get the methylation string for forward alignments.
    if (read_conversion == Constants.CONV_CT and direction == Constants.CONV_FORWARD) or \
     (read_conversion == Constants.CONV_GA and direction == Constants.CONV_REVERSE): #(CT and forward, CT_CT) or (GA and reversed, GA_CT)
        MI_string += ["M"] * 2
        bisulfite_sequence += ref_cdr
        for i in range(len_MI_string):
            #Get the identity of the first and second characters for determining the context.
            if MI_string[i + 1] != "M": 
                first_char = MI_string[i + 1]
            else:
                first_char = bisulfite_sequence[i + 1]
            if MI_string[i + 2] != "M":
                second_char = MI_string[i + 2]
            else:
                second_char = bisulfite_sequence[i + 2]
            #Call unmethylated cytosine
            if bisulfite_sequence[i] == 'T' and MI_string[i] == 'C':
                if first_char == 'G':
                    methylation_list.append('z')
                elif first_char == 'N' or second_char == 'N':
                    methylation_list.append('u')
                elif second_char == 'G':
                    methylation_list.append('x')
                else:
                    methylation_list.append('h')
            #Call methylated cytosine.
            elif bisulfite_sequence[i] == 'C' and MI_string[i] == 'M':
                if first_char == 'G':
                    methylation_list.append('Z')
                elif first_char == 'N' or second_char == 'N':
                    methylation_list.append('U')
                elif second_char == 'G':
                    methylation_list.append('X')
                else:
                    methylation_list.append('H')
            else:
                methylation_list.append('.')
    else:  #Get the methylation string for reverse alignments. (reversed and CT converted, CT_GA) OR (forward and GA converted, GA_GA)
        MI_string = MI_string[::-1] + ["M"] * 2
        bisulfite_sequence = bisulfite_sequence[::-1] + ref_cdr
        for i in range(len_MI_string):
            if MI_string[i + 1] != "M":
                first_char = MI_string[i + 1]
            else:
                first_char = bisulfite_sequence[i + 1]
            if MI_string[i + 2] != "M":
                second_char = MI_string[i + 2]
            else:
                second_char = bisulfite_sequence[i + 2]
            if bisulfite_sequence[i] == 'A' and MI_string[i] == 'G':
                if first_char == 'C':
                    methylation_list.append('z')
                elif first_char == 'N' or second_char == 'N':
                    methylation_list.append('u')
                elif second_char == 'C':
                    methylation_list.append('x')
                else:
                    methylation_list.append('h')
            elif bisulfite_sequence[i] == 'G' and MI_string[i] == 'M':
                if first_char == 'C':
                    methylation_list.append('Z')
                elif first_char == 'N' or second_char == 'N':
                    methylation_list.append('U')
                elif second_char == 'C':
                    methylation_list.append('X')
                else:
                    methylation_list.append('H')
            else:
                methylation_list.append('.')
        methylation_list = methylation_list[::-1]
    return "".join(methylation_list)


def writeSAMRecord(write_record, output_writer, **kwargs):
    """
    Given a write_record, the function write the SAM record to a SeqWriter
    in the output_writers dictionary.  The first position of the write_record 
    gives the correct SeqWriter to use.
    @param write_record: a tuple representing a record to be written the first
    position is the type of record, and the second position is the grouped SAM records
    @param output_writer: A writer that writes SAM records
    @param kwargs: named parameters that determine where to print to. (Not used)
    @return None
    """
    counter = 0
    for record in write_record[1]:
        counter += 1
        for subrecord in record:
            output_writer.write(subrecord)
    return counter
    

def performPostProcessFromMain(sam_files_list, pool, manager, outputfile, input_protocol, layout, 
                               command_line_string, postprocess_arguments):
    """
    The main entry point for BisPin post processing.  Takes the raw BFAST SAM files
    and replaces the converted reads with the bisulfite treated reads.  The MD tag 
    for the alignment is updated as well as other things.  May be done in parallel.
    @param sam_files_list: A list of strings that represent raw BFAST SAM files.
    @param pool: A pool of processes.  Possible set to None.  May not be used.
    @param manager: A manger server for multiprocessing.
    @param outputfile: A string representing where to write the output SAM file.
    @param input_protocol: A number from Constants indicating the input protocol.
    @param layout: A number from Constants indicating the type of layout.
    @param command_line_string: A string representing the command line parameters and arguments.  Added to the SAM file output.
    @param postprocess_arguments: A dictionary of extra arguments.
    @return A report string and the genome construction time.  
    """
    #Extract extra arguments.
    reads1 = postprocess_arguments.get(Constants.BISPIN_READS1, None)
    reads2 = postprocess_arguments.get(Constants.BISPIN_READS2, None)
    cutoff_bs = postprocess_arguments.get(Constants.BISPIN_CUTOFF_BS, Constants.BISPIN_DEFAULT_CUTOFF_BS)
    cutoff_r = postprocess_arguments.get(Constants.BISPIN_CUTOFF_R, Constants.BISPIN_DEFAULT_CUTOFF_R)
    gzip_switch = postprocess_arguments.get(Constants.GZIP, False)
    rescore_matrix = postprocess_arguments.get(Constants.BISPIN_RESCORE_MATRIX, None)
    start = postprocess_arguments.get(Constants.STARTREADNUM, None)
    end = postprocess_arguments.get(Constants.ENDREADNUM, None)
    single_process_switch = postprocess_arguments.get(Constants.SINGLEPROCESSPOST, True)
    genome_file = postprocess_arguments.get(Constants.FASTA, None)
    use_secondary_indexes = postprocess_arguments.get(Constants.USESECONDARYINDEXES, False)
    hairpin_post = postprocess_arguments.get(Constants.HAIRPINPOSTFILE, None)
    remove_comments = postprocess_arguments.get(Constants.HAIRPINPOSTFILE, False)
    #Set up either a double iterator or single iterator.
    if input_protocol == Constants.PROTOCOL_HAIRPIN or layout == Constants.LAYOUT_PAIRED:
        read_generator = SeqDoubleIterator.SeqDoubleIterator(reads1, reads2, file_type=Constants.FASTQ, gzip_switch = gzip_switch)
    else:
        read_generator = SeqIterator.SeqIterator(reads1, file_type=Constants.FASTQ, gzip_switch = gzip_switch)
    stats, genome_construction_time = selectBestAlignmentAndPostProcessSAMRecord(sam_files_list, 
                                                                                 pool, 
                                                                                 input_protocol, 
                                                                                 layout, 
                                                                                 read_generator, 
                                                                                 outputfile, 
                                                                                 cutoff_bs, 
                                                                                 cutoff_r, 
                                                                                 rescore_matrix, 
                                                                                 start, 
                                                                                 end, 
                                                                                 single_process_switch, 
                                                                                 gzip_switch, 
                                                                                 command_line_string, 
                                                                                 manager, 
                                                                                 genome_file,
                                                                                 use_secondary_indexes,
                                                                                 hairpin_post,
                                                                                 remove_comments,
                                                                                 )
    sys.stderr.write("%sWrote the SAM output...\n" % (logstr))
    sys.stderr.flush()
    report_string = makeOutputReport(stats)
    return report_string, genome_construction_time
    

def makeOutputReport(stats):
    """
    Create an output report summarizing alignment efficiency and methylation calling.
    @param stats: An object consisting of a list representing alignment efficiency statistics 
    and a dictionary representing methylation calling statistics.
    @return: A string representing an output report document.
    """
    alignment_breakdown = stats[0]
    uniquely_aligned_types = stats[1]
    methylation_stats = stats[2]
    total_alignments = (alignment_breakdown[0] + 0.0)
    if total_alignments == 0:
        output_string = ""
        output_string += "BisPin Alignment Report\n"
        output_string += "-------------------------\n"
        output_string += "No alignments were reported."
        return output_string
    output_string = ""
    output_string += "BisPin Alignment Report\n"
    output_string += "-------------------------\n"
    output_string += "Sequences analyzed in total: %s \n" % str(alignment_breakdown[0])
    output_string += "Number of alignments with a unique best hit:\t%s\t%s%%\n" % (str(alignment_breakdown[1]), 
                                                                                   str(100 * alignment_breakdown[1] / total_alignments)) 
    output_string += "Sequences that did not map uniquely:\t%s\t%s%%\n" % (str(alignment_breakdown[2]), 
                                                                           str(100 * alignment_breakdown[2] / total_alignments)) 
    try:
        some_percent = str(100 * alignment_breakdown[5] / (0.0 + alignment_breakdown[5] + alignment_breakdown[2]))
    except ZeroDivisionError:
        some_percent = str(None)
    output_string += "Sequences that mapped uniquely after rescoring:\t%s\t%s%%\t%s%%\n" % (str(alignment_breakdown[5]), 
                                                                                            str(100 * alignment_breakdown[5] / total_alignments), some_percent)
    output_string += "Sequences that did not map because no hits were found:\t%s\t%s%%\n" % (str(alignment_breakdown[3]), 
                                                                                             str(100 * alignment_breakdown[3] / total_alignments))
    output_string += "Sequences that did not map because they were filtered out:\t%s\t%s%%\n" % (str(alignment_breakdown[4]), 
                                                                                                 str(100 * alignment_breakdown[4] / total_alignments))
    output_string += "\n"
    output_string += "Breakdown of uniquely aligned sequences.\n"
    uat_keys = list(uniquely_aligned_types.keys())
    uat_keys.sort()
    for key in uat_keys:
        output_string += "%s:\t%s\t%s\n" % (key, str(uniquely_aligned_types[key]), 
                                            str(Constants.CONV_DESC.get(key, "DESCRIPTION_NOT_FOUND")))
    output_string += "\n"
    output_string += calculateMethylationStats(methylation_stats)
    return output_string
    
    
def calculateMethylationStats(methylation_stats):
    """
    Produces a string document summarizing the methylation calling.
    @param methylation_stats: A dictionary counting each methylation calling context found in the post processed SAM alignments.
    @return A string report summarizing the methylation calling.
    """
    output_string = ""
    Z = methylation_stats.get('Z', 0) #Z - methylated C in the CpG context
    z = methylation_stats.get('z', 0) #z - unmethylated C in the CpG context
    X = methylation_stats.get('X', 0) #X - methylated C in the CHG context
    x = methylation_stats.get('x', 0) #x - unmethylated C in the CHG context
    H = methylation_stats.get('H', 0) #H - methylated C in the CHH context
    h = methylation_stats.get('h', 0) #h - unmethylated C in the CHH context
    U = methylation_stats.get('U', 0)
    u = methylation_stats.get('u', 0)
    totalCs = Z + z + X + x + H + h + u + U
    if z + Z != 0:
        CpG = 100 * Z / (z + Z + 0.0)
        CpG_p = 100 * z / (z + Z + 0.0)
    else:
        CpG = None
        CpG_p = None
    if x + X != 0:
        CHG = 100 * X / (x + X + 0.0)
        CHG_p = 100 * x / (x + X + 0.0)
    else:
        CHG = None
        CHG_p = None
    if h + H != 0: 
        CHH = 100 * H / (h + H + 0.0)
        CHH_p = 100 * h / (h + H + 0.0)
    else:
        CHH = None
        CHH_p = None
    if u + U != 0:
        CN = 100 * U / (u + U + 0.0)
        CN_p = 100 * u / (u + U + 0.0)
    else:
        CN = None
        CN_p = None
    methylation_stats['totalCs'] = totalCs
    methylation_stats['CpG'] = CpG
    methylation_stats['CHG'] = CHG
    methylation_stats['CHH'] = CHH
    methylation_stats['CN'] = CN
    output_string += "Final Cytosine Methylation Report\n"
    output_string += "----------------------------------\n"
    output_string += "Total number of C's analyzed:\t%s\n" % (str(totalCs))
    output_string += "\n"
    output_string += "Total methylated C's in the CpG context:\t%s\t%s%%\n" % (str(Z), str(CpG))
    output_string += "Total methylated C's in the CHG context:\t%s\t%s%%\n" % (str(X), str(CHG))
    output_string += "Total methylated C's in the CHH context:\t%s\t%s%%\n" % (str(H), str(CHH))
    output_string += "Total methylated C's in the CNN context:\t%s\t%s%%\n" % (str(U), str(CN))
    output_string += "\n"
    output_string += "Total C to T conversions in the CpG context:\t%s\t%s%%\n" % (str(z), str(CpG_p))
    output_string += "Total C to T conversion in the CHG context:\t%s\t%s%%\n" % (str(x), str(CHG_p))
    output_string += "Total C to T conversion in the CHH context:\t%s\t%s%%\n" % (str(h), str(CHH_p))
    output_string += "Total C to T conversion in the CNN context:\t%s\t%s%%\n" % (str(u), str(CN_p))
    return output_string


def createRescoreMatrix(score_lines, methylMatrix = False):
    """
    Creates the rescoring matrix function.
    @param score_lines: A string representing the rescoring matrix function.
    @return: A dictionary representing the rescoring matrix function.
    """
    score_matrix = {}
    score_lines = [x for x in score_lines if x != '']
    columns = score_lines[0].upper().strip().split()
    open_scores = score_lines[6].strip().split(',')
    if len(open_scores) == 1:
        open_scores *= 2
    extension_scores = score_lines[7].strip().split(',')
    if len(extension_scores) == 1:
        extension_scores *= 2
    score_matrix[Constants.OPEN_DEL] = float(open_scores[0])
    score_matrix[Constants.OPEN_INS] = float(open_scores[1])
    score_matrix[Constants.EXT_DEL] = float(open_scores[0])
    score_matrix[Constants.EXT_INS] = float(open_scores[1])
    for i in range(1, 6):
        row = score_lines[i].strip().upper().split()
        row_dict = {}
        for i in range(1, len(row)):
            row_dict[columns[i - 1]] = float(row[i])
        score_matrix[row[0]] = row_dict
    if not methylMatrix:
        methylation_call_penalty = (score_matrix['C']['C'] + score_matrix['A']['A'] + 
                                    score_matrix['G']['G'] + score_matrix['T']['T']) / 4.0
        score_matrix["GA"] = methylation_call_penalty
        score_matrix["CT"] = methylation_call_penalty
    else:
        score_matrix["GA"] = score_matrix["G"]["A"]
        score_matrix["CT"] = score_matrix["C"]["T"]
    return score_matrix
        

def handleRescoreOptions(rescore, rescore_matrix_1, rescore_matrix_2):
    """
    Returns the rescoring matrix dictionary.  Parses default settings.
    @param rescore: A boolean determining if the rescoring will be done.
    @param rescore_matrix_1: A string indicating the location of a file representing the rescoring matrix.
    @return: The rescoring matrix dictionary or None if rescore is False. 
    """
    if rescore and rescore_matrix_1 != None:
        if rescore_matrix_2 == None:
            rescore_matrix_1 = createRescoreMatrix(open(rescore_matrix_1).readlines())
        else:
            rescore_matrix_1 = createRescoreMatrix(open(rescore_matrix_1).readlines(), methylMatrix = True)
            rescore_matrix_2 = createRescoreMatrix(open(rescore_matrix_2).readlines(), methylMatrix = True)
    elif rescore and rescore_matrix_1 == None:
        rescore_matrix_1 = createRescoreMatrix(Constants.RESCORE_DEFAULT.split("\n"))
    else:
        rescore_matrix_1 = None
    return (rescore_matrix_1, rescore_matrix_2)


def main():
    """
    The main function is the entry point for the module.  This function creates the command line interface and parses the input.  It displatches the postprocessing functions.
    """
    now = datetime.datetime.now()
    usage = "usage: %prog [options] <reference_genome_file> <sam_input_directory> <output_file> { <reads1> <reads2> | <reads> }"
    version = "%prog " + Constants.version
    description = "This program postprocesses SAM files for bisulfite treatment.  The SAM files in the <sam_input_directory> need to have the string CT_CT or CT_GA or GA_CT or GA_GA in their name for correct processing, and they must end '.sam'.  The <output_file> is the SAM file where the results will be written to.  It will be used as a prefix to write an alignment report file.  Two reads files are specified for hairpin and paired end data in <reads1> and <reads2>, but single end data uses one file <reads>.  The T-enriched file is the first file and the A-enriched file is the second file."
    epilog = Constants.creation_string
    p = optparse.OptionParser(usage = usage, 
                              version = version, 
                              description = description, 
                              epilog = epilog, 
                              formatter = IndentedHelpFormatterWithNL.IndentedHelpFormatterWithNL())
    p.add_option('--protocol', '-p', help='This gives the construction protocol that will be used. [default: %default]\n' + 
                 str(Constants.PROTOCOL_DIRECTIONAL) + ': directional\n' + 
                 str(Constants.PROTOCOL_BIDIRECTIONAL) + ': bidirectional\n' + 
                 str(Constants.PROTOCOL_PBAT) + ': PBAT\n' +
                 str(Constants.PROTOCOL_HAIRPIN) + ': hairpin\n' ,
                 #str(Constants.PROTOCOL_ONLYCT) + ': only CT reads against a CT genome\n' +
                 #str(Constants.PROTOCOL_ONLYGA) + ': only GA reads against a GA genome\n' ,
                 default=Constants.PROTOCOL_DIRECTIONAL)
    p.add_option('--gzip', '-z', help='The output file will be gzip compressed. [default: %default]', action='store_true', default=False)
    p.add_option('--startReadNum', '-s', help='Specifies the read to begin with.  Reads are numbered starting at 0. (skip the first startReadNum-1 reads) [default: %default]', default = None)
    p.add_option('--endReadNum', '-e', help='Specifies the last read to use (inclusive) [default: %default]', default = None)
    p.add_option('--removeComments', '-C', help='Do not include comments at the top of the SAM file.  This can be used while partitioning the input so that the output SAM files can be easily concatenated.', action='store_true', default=False)
    p.add_option('--multiProcessPost', '-W', help='Use multiple processes for post processing.  This multiprocessing uses six processes, so this switch should only be engaged if the machine has six or more cores.', action = 'store_true', default = False)
    p.add_option('--useSecondaryIndexes', '-i', help = "Turn this flag on if the SAM files were created with secondary indexes." , action='store_true', default=False)
    filter_group = optparse.OptionGroup(p, "Filter Options", "Specify the function for filtering out low quality alignments and reporting the reads as unmapped.")
    filter_group.add_option('--cutoff_bs', '-b', help='Function for filtering low quality alignments for bisulfite treated reads. [default: %default]', default=Constants.BISPIN_DEFAULT_CUTOFF_BS)
    filter_group.add_option('--cutoff_r', '-r', help='Function for filtering low quality alignments for recovered reads.  This is for hairpin data only. [default: %default]', default=Constants.BISPIN_DEFAULT_CUTOFF_R)
    p.add_option_group(filter_group)
    rescore_group = optparse.OptionGroup(p, "Rescore Options", "Ambiguously aligned reads can be rescored to see if one will be unique.")
    rescore_group.add_option('--noRescore', '-N', help='Do NOT rescore ambiguously aligned reads to find unique alignments. [default: %default]', action='store_true', default = False)
    rescore_group.add_option('--rescore_matrix_1', '-x', help='The location of the rescoring matrix used to rescore ambiguously aligned reads for a unique score.  If the second rescoring matrix is specified, it is assumed that this rescoring matrix will be for C to T methylation rescoring. The default is:\n' + Constants.RESCORE_DEFAULT, default=None)
    rescore_group.add_option('--rescore_matrix_2', '-X', help='The location of the rescoring matrix used to rescore ambiguously aligned reads for G to A methylation rescoring for a unique score.  [default: %default]', default=None)
    p.add_option_group(rescore_group)
    options, args = p.parse_args()
    if len(args) == 0:
        p.print_help()
    if len(args) < 4 or len(args) > 5:
        p.error('Not enough arguments.  Check the usage.')
    try:
        input_protocol = int(options.protocol)
    except ValueError:
        p.error("The protocol is not a number.")
    if input_protocol < Constants.PROTOCOL_DIRECTIONAL or input_protocol > Constants.PROTOCOL_ONLYGA:
        p.error('The protocol number is invalid.  Please check.')
    sys.stderr.write("%sStarting BisPin_postprocess.  Current time: %s \n" % (logstr, str(now)))
    rescore_matrix = handleRescoreOptions(not options.noRescore, options.rescore_matrix_1, options.rescore_matrix_2)
    command_line_string = "options: %s args: %s rescore_matrix: %s" % (str(options), str(args), str(rescore_matrix))
    sys.stderr.write("%s Command line -- %s\n" % (logstr, command_line_string))
    cutoff_bs = float(options.cutoff_bs)
    cutoff_r = float(options.cutoff_r)
    layout = Constants.LAYOUT_SINGLE
    reads2 = None
    if len(args) == 5:
        reads2 = args[4]
        layout = Constants.LAYOUT_PAIRED
    start = None if options.startReadNum == None else int(options.startReadNum)
    end = None if options.endReadNum == None else int(options.endReadNum)
    if (start != None or end != None) and input_protocol == Constants.PROTOCOL_HAIRPIN and options.useSecondaryIndexes == True:
        p.error("The partitioning feature for the hairpin protocol is not supported at this time for secondary indexes.  The Linux head and tail commands can be used to partition data.")
    postprocess_arguments = {Constants.FASTA: args[0],
                             Constants.BISPIN_READS1: args[3], 
                             Constants.BISPIN_READS2: reads2,
                             Constants.BISPIN_CUTOFF_BS: cutoff_bs, 
                             Constants.BISPIN_CUTOFF_R: cutoff_r,
                             Constants.BISPIN_RESCORE_MATRIX: rescore_matrix,
                             Constants.GZIP: options.gzip,
                             Constants.STARTREADNUM: start,
                             Constants.ENDREADNUM: end,
                             Constants.SINGLEPROCESSPOST : not options.multiProcessPost,
                             Constants.USESECONDARYINDEXES : options.useSecondaryIndexes,
                             Constants.REMOVECOMMENTS : options.removeComments,
                             }
    outputfile = args[2]
    sam_files_dir = args[1]
    if not os.path.isdir(sam_files_dir):
        p.error("The sam directory is invalid.  Please ensure that it is a directory.")
    onlyfiles = [f for f in os.listdir(sam_files_dir) if os.path.isfile(os.path.join(sam_files_dir, f))]
    onlyfiles = [os.path.join(sam_files_dir, f) for f in onlyfiles]
    onlyfiles = [f for f in onlyfiles if Constants.SAM.lower() in os.path.splitext(f)[-1].lower()]
    pool = None
    manager = multiprocessing.Manager()
    report_string, genome_counstruction_time = performPostProcessFromMain(onlyfiles, pool, manager, 
                                                                          outputfile, input_protocol, 
                                                                          layout, command_line_string, 
                                                                          postprocess_arguments)
    later = datetime.datetime.now()
    timing_report = BisPin_util.createTimingInfoReport(None, None, None, None, None, now, later, genome_counstruction_time)
    final_report = report_string + "\n" + timing_report + "\n\n" + command_line_string
    open(os.path.join(outputfile + '.BisPin.report'), 'w').write(final_report)
    sys.stdout.write(logstr + final_report + "\n")
    sys.stderr.write("%sEnding BisPin_postprocess.  Current time: %s Elapsed time: %s \n" % (logstr, str(later), str(later - now)))


if __name__ == "__main__":
    main()