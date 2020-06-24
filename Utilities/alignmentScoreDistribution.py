#!/usr/bin/python
import datetime
import optparse
import sys
from . import IndentedHelpFormatterWithNL 
from . import Constants
from . import SeqIterator
import math

"""
This file calculates a list of alignment scores for a SAM file.
The user can either specify a number of buckets for a histogram or have the results returned as a list of scores.
The final score is given as the SAM file alignment score divided by the sequence length of the SAM record.
Useful statistics including the max, the min, the average, and the median are given.
@author: Jacob Porter
"""

bucket_count_default = -1

def isUnmapped(flag):
    """Checks if the SAM flag indicates an unmapped read."""
    return ((int(flag) >> 2) % 2) == 1

def isFiltered(flag):
    """Checks if the SAM flag indicates a filtered reads."""
    return ((int(flag) >> 9) % 2) == 1

def findASDistribution(sam_file, use_edit_distance, bucket_count = bucket_count_default):
    sam_input = SeqIterator.SeqIterator(sam_file, file_type=Constants.SAM)
    as_max = -sys.maxsize - 1
    as_min = sys.maxsize
    no_hits = 0
    total_records = 0
    scores = []
    for record in sam_input:
        if record[Constants.SAM_KEY_RNAME].startswith(Constants.SAM_VALUE_STAR) or isUnmapped(record[Constants.SAM_KEY_FLAG]) :
            no_hits += 1
            continue
        if not use_edit_distance:
            alignment_score = float(record[Constants.SAM_KEY_ALIGNMENT_SCORE]) / float(len(record[Constants.SAM_KEY_SEQ]))
        else:
            alignment_score = float(record[Constants.SAM_KEY_DISTANCE]) / float(len(record[Constants.SAM_KEY_SEQ]))
        scores.append(alignment_score)
        if alignment_score > as_max:
            as_max = alignment_score
        if alignment_score < as_min:
            as_min = alignment_score
        total_records += 1
    scores.sort()
    total_score = sum(scores)
    avg_score = total_score / (total_records + 0.0)
    median_position = int(total_records / 2)
    median_score = scores[median_position]
    if bucket_count > 0:
        buckets = [0] * (bucket_count + 1)
        as_range = as_max - as_min
        bucket_length = as_range / bucket_count
        score_bounds = as_min
        score_bounds_list = []
        for _ in range(bucket_count):
            score_bounds += bucket_length
            score_bounds_list.append(score_bounds)
        for alignment_score in scores:
            bucket_index = int(math.floor((alignment_score - as_min) / bucket_length))
            buckets[bucket_index] += 1
        freq = list(map((lambda x: x / (total_records + 0.0)), buckets))
        return (freq, score_bounds_list, total_records, as_max, as_min, buckets, avg_score, median_score, no_hits)
    else:
        return (scores, "", total_records, as_max, as_min, "", avg_score, median_score, no_hits)

def main():
    """
    Calculates a distribution of the alignment scores for a SAM file.  The distribution is printed to standard out.
    """
    now = datetime.datetime.now()
    usage = "usage: %prog [options] <sam_file> "
    version = "%prog " + Constants.version
    description = ""
    epilog = ""
    p = optparse.OptionParser(usage = usage, version = version, 
                              description = description, epilog = epilog, 
                              formatter = IndentedHelpFormatterWithNL.IndentedHelpFormatterWithNL())
    p.add_option('--buckets','-b',help='The number of buckets to produce in the frequency distribution. If this is set to zero or a negative value, then a list of scores will be returned.  [default: %default]', default = bucket_count_default)
    p.add_option('--edit','-e',help='Use the edit distance instead of the alignment score.  [default: %default]', action='store_true', default = False)
    p.add_option()
    options, args = p.parse_args()
    if len(args) != 1:
        p.print_help()
        p.error("There must be one argument.")
    sam_file = args[0]
    sys.stdout.write("Processing the file %s with %s buckets.\n" % (sam_file, str(options.buckets)))
    stats = findASDistribution(sam_file, options.edit, bucket_count = int(options.buckets))
    (freq, score_bounds_list, total_records, as_max, as_min, buckets, avg_score, median_score, no_hits) = stats
    stats = [str(x) for x in stats]
    result_string = "Frequency distribution / scores: %s\nScore bounds: %s\nReads processed: %s\nMaximum score: %s\nMinimum score: %s\nCounts: %s\nAverage score: %s\nMedian score: %s\nNumber of no hit alignments: %s\n" % (freq, score_bounds_list, total_records, as_max, as_min, buckets, avg_score, median_score, no_hits)
    later = datetime.datetime.now()
    runtime = later - now
    sys.stdout.write(result_string + "Runtime: " + str(runtime) + "\n")

if __name__ == "__main__":
    main()