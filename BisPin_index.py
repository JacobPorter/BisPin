#!/usr/bin/python
import optparse, os, sys, multiprocessing, subprocess, datetime
from Utilities import Constants
from Utilities import IndentedHelpFormatterWithNL
from Utilities import BisPin_util

logstr = "\nBisPin_index: "

"""
@author: Jacob Porter
@summary: Builds BFAST indices.  The reference genomes must already be converted.
TODO: Add depth, etc. to indices.
TODO: Creating the brg file should probably be separated because it needs to be done before the index file is created.
"""

def fasta2brg(path_to_bfast, directory, filename, CTorGA, orig = False):
    """
    Converts a FASTA file into a BFAST binary format, brg.
    This function is called after the reference genome is converted.
    @param path_to_bfast: a string giving the location of the BFAST executable.
    @param directory: the directory where the reference FASTA file is located
    @param filename: the filename of the reference FASTA file
    @param CTorGA: if True, indicates that the FASTA file is C to T converted
    If False, the FASTA file is G to A converted.  
    This is used to get the correct name of the reference file.
    @param orig: when True, the FASTA file is an unconverted reference file.
    """
    if not orig:
        str1 = Constants.CONV_CT if CTorGA else Constants.CONV_GA
        genome_path = os.path.join(directory, filename + ".BisPin." + str1)
    else:
        genome_path = os.path.join(directory, filename)
    BRG_path = os.path.join(genome_path + ".nt" + ".brg")
    if os.path.exists(BRG_path):
        sys.stderr.write("%sThe BRG file %s already exists, so it will not be created.\n" % (logstr, str(BRG_path)))
        sys.stderr.flush()
    else:
        sys.stderr.write("%sThe BRG file %s does NOT exist, so it will be created.\n" % (logstr, str(BRG_path)))
        sys.stderr.flush()
        subprocess.call([path_to_bfast, 'fasta2brg', '-f', genome_path, "-A", "0"])

def create_index(fastafile, path_to_bfast, mask, hashWidth, indexNumber, depth, num_threads, tmpDir, CTorGA):
    """
    Create a BFAST index
    @see: BFAST documentation for more information on the index building paraemters.
    @param fastafile: A string giving the location of a FASTA file. 
    @param path_to_bfast: a string giving the location of the BFAST executable.
    @param mask: A BFAST index mask
    @param hashWidth: A BFAST index hash width
    @param indexNumber: A BFAST index number
    @param depth: A BFAST depth
    @param num_threads: the number of threads for BFAST to use (This seems to have no effect.)
    @param tmpDir: A temporary directory for BFAST to write to.
    @param CTorGA: When True, the FASTA file is expected to conform to the naming conventions of a 
    C to T converted reference genome.  When False, it should correspond to a G to A converted file.
    """
    str1 = "CT" if CTorGA else "GA"
    converted_genome_path = os.path.join(fastafile + ".BisPin." + str1)
    converted_genome_path = str(converted_genome_path)
    subprocess.call([path_to_bfast, 'index', '-f', converted_genome_path, '-m', mask, '-w', hashWidth, '-i', indexNumber, '-d', depth, '-n', num_threads, '-T', tmpDir])

def create_regular_index(fastafile, path_to_bfast, mask, hashWidth, indexNumber, depth, num_threads, tmpDir):
    """
    Creates a BFAST index for a regular unconverted reference genome.
    @see: create_index help documentation has information on the parameters.
    """
    genome_path = os.path.join(fastafile)
    subprocess.call([path_to_bfast, 'index', '-f', genome_path, "-m", mask, "-w", hashWidth, '-i', indexNumber, '-d', depth, '-n', num_threads, '-T', tmpDir])

def create_bwa_converted_index(fastafile, path_to_bwa, block_size, CTorGA):
    """
    Experimental code for BWA indices.
    @status: Experimental 
    """
    str1 = "CT" if CTorGA else "GA"
    converted_genome_path = os.path.join(fastafile + ".BisPin." + str1)
    subprocess.call([path_to_bwa, 'index', '-b', block_size, converted_genome_path])

def create_bwa_regular_index(fastafile, path_to_bwa, block_size):
    """
    Experimental code for BWA indices.
    @status: Experimental
    """
    genome_path = os.path.join(fastafile)
    subprocess.call([path_to_bwa, 'index', '-b', block_size, genome_path])

def main():
    """
    The entry point into the program.
    """
    #Create the command line interface.
    usage = "usage: %prog [options] <reference_genome_file> "
    version = "%prog " + Constants.version
    description = "Create indices for the reference genome.  The file naming convention for the converted genome is expected to be in BisPin_covert format."
    epilog = Constants.creation_string
    p = optparse.OptionParser(usage = usage, version = version, description = description, 
                              epilog = epilog, formatter = IndentedHelpFormatterWithNL.IndentedHelpFormatterWithNL())
    p.add_option('--path', '-P', help='The path to the BFAST executable file.  If this option is not used, then the PATH variable will be searched for the executable.', default = None)
    p.add_option('--tmpDir','-T',help= 'Specifies the directory to store temporary files.  [default: the directory where the outputfile is located.]', default = None)
    p.add_option('--regular', '-r', help='Create an index for the regular unconverted reference genome. This is for processing hairpin data for the recovery step.  [default: %default]', action='store_true', default=False)
    #p.add_option('--aligner', '-a', help='The code string for the aligner to use.  [default: %default]\nOptions:\n' + Constants.BFAST +'\n' + Constants.BWA, default=Constants.BFAST)
    p.add_option('--sequential', '-s', help='Create the indices for the converted reference genome in sequence instead of in parallel. [default: %default]', action='store_true', default=False)
    bfast_group = optparse.OptionGroup(p, "BFAST index options.")
    bfast_group.add_option('--mask', '-m', help='The mask or spaced seed to use. The mask is a set of zero and ones (must start and end with a one). [default: %default]', default = '11111111111111111111')
    bfast_group.add_option('--hashWidth', '-w', help='The length of the hashed string (key size for the hash index). This must be less than or equal to the number of ones in the mask. [default: %default]', default = 14)
    bfast_group.add_option('--indexNumber', '-i', help='Specifies this is the ith index you are creating. This is useful when multiple indexes from the same reference are to be created (in the same space). [default: %default]', default='1')
    #p.add_option('--num_threads', '-n', help='The number of threads for each BFAST process to use. [default: %default]', default='1')
    #p.add_option('--depth' , '-d', help= 'The depth of splitting (d).  The index will be split into 4^d parts. [default: %default]', default='0')
    p.add_option_group(bfast_group)
#     bwa_group = optparse.OptionGroup(p, "BWA index options.")
#     p.add_option('--block_size', '-b', help='block size for the bwtsw algorithm [default: %default]', default='10000000')
#     p.add_option_group(bwa_group)
    #Extract arguments
    options, args = p.parse_args()
    now = datetime.datetime.now()
    if len(args) == 0:
        p.print_help()
    if len(args) != 1:
        p.error("There must be one argument.")
    fastafile = args[0]
    #path_to_aligner = args[0]
    aligner_type = Constants.BFAST #options.aligner.upper()
    if aligner_type != Constants.BFAST and aligner_type != Constants.BWA:
        p.error("The aligner is not recognized.  Please choose the aligner from the given options.")
    if options.path == None:
        path_to_aligner = BisPin_util.which("bfast")
    else:
        path_to_aligner = BisPin_util.which(options.path)
    if path_to_aligner == None:
        p.error("The BFAST executable could not be found.  Please check the path.")
    path_to_aligner = str(path_to_aligner)
    try:
        hashWidth = int(options.hashWidth)
    except ValueError:
        p.error("The hash width argument must be a positive integer.")
    regular = options.regular
    if not os.path.exists(fastafile):
        p.error("The file at %s could not be located." % fastafile)
    if not os.path.exists(path_to_aligner):
        p.error("The aligner could not be found at %s" % path_to_aligner)
    mask = options.mask
    hashWidth = str(hashWidth)
    indexNumber = str(options.indexNumber)
    sequential = options.sequential
    tmpDir = options.tmpDir
    num_threads = '1'#options.num_threads  # Multithreading for constructing indices does not seem to do anyting.
    depth = '0'#options.depth #Unsupported right now.
    filename = os.path.basename(fastafile)
    directory = os.path.dirname(fastafile)
    if tmpDir == None:
        tmpDir = directory
    tmpDir = str(tmpDir)
    if not tmpDir.endswith("/"):
        tmpDir += "/"
    #Create indices
    sys.stderr.write("%sCreating indices for %s\n" % (logstr, fastafile))
    sys.stderr.flush()
    if aligner_type == Constants.BFAST: #Create BFAST indices
        sys.stderr.write(logstr + "Converting the FASTA format into the BFAST binary version (BRG)...\n")
        sys.stderr.flush()
        finish_string = "%sFinished creating %s index with mask %s and hash width %s for %s with\n" + aligner_type
        if not regular:
            p_CT_bin = multiprocessing.Process(target=fasta2brg, args=(path_to_aligner, directory, filename, True))
            p_CT_bin.start()
            p_GA_bin = multiprocessing.Process(target=fasta2brg, args=(path_to_aligner, directory, filename, False))
            p_GA_bin.start()
            p_CT_bin.join()
            p_GA_bin.join()
            sys.stderr.write(logstr + "Finished creating the BFAST BRG files!\n")
            p_CT = multiprocessing.Process(target=create_index, args=(fastafile, path_to_aligner, mask, hashWidth, indexNumber, depth, num_threads, tmpDir, True))
            p_CT.start()
            if sequential:
                p_CT.join()
                sys.stderr.write(finish_string % (logstr, "C to T", mask, hashWidth, fastafile))
            p_GA = multiprocessing.Process(target=create_index, args=(fastafile, path_to_aligner, mask, hashWidth, indexNumber, depth, num_threads, tmpDir, False))
            p_GA.start()
            if not sequential:
                p_CT.join()
                sys.stderr.write(finish_string % (logstr, "C to T", mask, hashWidth, fastafile))
            p_GA.join()
            sys.stderr.write(finish_string % (logstr, "G to A", mask, hashWidth, fastafile))
        else:
            p_orig_bin = multiprocessing.Process(target=fasta2brg, args=(path_to_aligner, directory, filename, True, True))
            p_orig_bin.start()
            p_orig_bin.join()
            sys.stderr.write(logstr + "Finished creating the BFAST BRG files!\n")
            p_orig = multiprocessing.Process(target=create_regular_index, args=(fastafile, path_to_aligner, mask, hashWidth, indexNumber, depth, num_threads, tmpDir))
            p_orig.start()
            p_orig.join()
            sys.stderr.write(finish_string % (logstr, "original", mask, hashWidth, fastafile))
    elif aligner_type == Constants.BWA: #Create BWA indices
        finish_string = "%sFinished creating the %s index with " + aligner_type 
        if not regular:
            p_CT = multiprocessing.Process()
            p_CT.start()
            if sequential:
                p_CT.join()
                sys.stderr.write(finish_string % (logstr, "C to T"))
            p_GA = multiprocessing.Process()
            p_GA.start()
            if not sequential:
                p_CT.join()
                sys.stderr.write(finish_string % (logstr, "C to T"))
            p_GA.join()
            sys.stderr.write(finish_string % (logstr, "G to A"))
        else:
            p_orig = multiprocessing.Process()
            p_orig.start()
            p_orig.join()
            sys.stderr.write(finish_string % (logstr, "original"))
    sys.stderr.flush()
    later = datetime.datetime.now()
    sys.stderr.write("%sElapsed time -- %s\n" % (logstr, str(later - now)))
    
    
if __name__ == '__main__':
    main()
