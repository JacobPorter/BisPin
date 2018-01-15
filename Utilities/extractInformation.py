#!/usr/bin/python
"""
Utility files for extracting information for accuracy files and from
mapper program result files.
@author: Jacob Porter
"""
import re
import sys

def splitAccuracyLine(line):
    try:
        value = float(line.strip().split(":")[1])
    except ValueError:
        value = -1
    except TypeError:
        value = -1
    return value

def extractAccuracy(file_location):
    fd = open(file_location, 'r')
    value_dictionary = {}
    value_dictionary["Total"] = 0.0
    value_dictionary["Recall"]= 0.0
    value_dictionary["Precision"] = 0.0
    value_dictionary["F1_Score"] = "NaN"
    for line in fd:
        if line.startswith("Total number of reads from the input"):
            value_dictionary["Total"] = splitAccuracyLine(line)
        elif line.startswith("Percent correct of total"):
            value_dictionary["Recall"] = splitAccuracyLine(line)
        elif line.startswith("Percent correct of analyzed"):
            value_dictionary["Precision"] = splitAccuracyLine(line)
    if value_dictionary["Recall"] != 0.0 or value_dictionary["Precision"] != 0.0:
        value_dictionary["F1_Score"] = 2 * value_dictionary["Recall"] * value_dictionary["Precision"] / (value_dictionary["Precision"] + value_dictionary["Recall"])
    return value_dictionary

def splitBisPinLine(line):
    line = line.replace("%", "")
    line_list = line.strip().split(":")
    line_list = line_list[1].split()
    #sys.stderr.write("BisPin line_list: %s\n" % str(line_list))
    try:
        second = float(line_list[1])
        second /= 100
    except IndexError:
        second = 1.0
    except ValueError:
        second = None
    return (int(line_list[0]), second)

def extractBisPin(file_location):
    fd = open(file_location, 'r')
    value_dictionary = {}
    for line in fd:
        if line.startswith("Sequences analyzed in total"):
            value_dictionary["Reads"] = splitBisPinLine(line)[0]
        elif line.startswith("Number of alignments with a unique best hit"):
            value_dictionary["Unique"] = splitBisPinLine(line)[1]
        elif line.startswith("Sequences that did not map uniquely"):
            value_dictionary["Ambig"] = splitBisPinLine(line)[1]
        elif line.startswith("Sequences that did not map because no hits were found"):
            value_dictionary["Unmap"] = splitBisPinLine(line)[1]
        elif line.startswith("Sequences that did not map because they were filtered out"):
            value_dictionary["Filt"] = splitBisPinLine(line)[1]
        elif "methylated C's in the CpG" in line:
            value_dictionary["CpG"] = splitBisPinLine(line)[1]
        elif "methylated C's in the CHG" in line:
            value_dictionary["CHG"] = splitBisPinLine(line)[1]
        elif "methylated C's in the CHH" in line:
            value_dictionary["CHH"] = splitBisPinLine(line)[1]
    value_dictionary["Unmap + Filt"] = value_dictionary["Unmap"] + value_dictionary["Filt"]
    return value_dictionary

def splitBismarkLine(line, total):
    line_list = line.strip().split(":")
    if total == None:
        return int(line_list[1])
    elif total < 0:
        return (None, float(line_list[1].replace("%", "")) / 100.0)
    else:
        return (int(line_list[1]), int(line_list[1]) / (total + 0.0))

def extractBismark(file_location):
    fd = open(file_location, 'r')
    value_dictionary = {}
    for line in fd:
        if line.startswith("Sequences analysed in total"):
            total_reads = splitBismarkLine(line, None)
            value_dictionary["Reads"] = (total_reads, 1.0)[0]
            break
    fd.close()
    fd = open(file_location, 'r')
    for line in fd:
        if line.startswith("Number of alignments with a unique"):
            value_dictionary["Unique"] = splitBismarkLine(line, total_reads)[1]
        elif line.startswith("Sequences did not map"):
            value_dictionary["Ambig"] = splitBismarkLine(line, total_reads)[1]
        elif line.startswith("Sequences with no alignments"):
            value_dictionary["Unmap + Filt"] = splitBismarkLine(line, total_reads)[1]
        elif "methylated in CpG" in line:
            value_dictionary["CpG"] = splitBismarkLine(line, -1)[1]
        elif "methylated in CHG" in line:
            value_dictionary["CHG"] = splitBismarkLine(line, -1)[1]
        elif "methylated in CHH" in line:
            value_dictionary["CHH"] = splitBismarkLine(line, -1)[1]
    return value_dictionary

def splitBWAMethLine(line):
    line = line.strip()
    line_list = re.split(":|\t", line)
    line_list = [item for item in line_list if item != '']
    #sys.stderr.write("BWA line_list: %s\n" % str(line_list))
    try:
        second = float(line_list[2])
    except IndexError:
        second = 1.0
    except ValueError:
        second = float('nan')
    return (int(line_list[1]), second)

def extractBWAMeth(file_location):
    """
    Reads:  491106
    Unique: 431513  0.878656
    Ambig:  48559   0.098877
    Unmap:  11030   0.022460
    Filt:   4       0.000008
    Unmap + Filt:   11034   0.022468
    """
    fd = open(file_location, 'r')
    value_dictionary = {}
    for line in fd:
        if line.startswith("Reads"):
            value_dictionary["Reads"] = splitBWAMethLine(line)[0]
        elif line.startswith("Unique"):
            value_dictionary["Unique"] = splitBWAMethLine(line)[1]
        elif line.startswith("Ambig"):
            value_dictionary["Ambig"] = splitBWAMethLine(line)[1]
        elif line.startswith("Unmap:"):
            value_dictionary["Unmap"] = splitBWAMethLine(line)[1]
        elif line.startswith("Filt:"):
            value_dictionary["Filt"] = splitBWAMethLine(line)[1]
        elif line.startswith("Unmap + Filt"):
            value_dictionary["Unmap + Filt"] = splitBWAMethLine(line)[1]
    return value_dictionary

def splitWaltLine(line):
    line = line.replace("[", "").replace("]", "").replace("(", ":").replace(")", "").replace("%", "")
    line_list = line.strip().split(":")
    try:
        second = float(line_list[2])
    except IndexError:
        second = 100.0
    second /= 100.0
    return (int(line_list[1]), second)

def extractWalt(file_location):
    fd = open(file_location, 'r')
    value_dictionary = {}
    for line in fd:
        if "TOTAL" in line:
            value_dictionary["Reads"] = splitWaltLine(line)[0]
        elif "UNIQUELY" in line:
            value_dictionary["Unique"] = splitWaltLine(line)[1]
        elif "AMBIGUOUS" in line:
            value_dictionary["Ambig"] = splitWaltLine(line)[1]
        elif "UNMAPPED" in line:
            value_dictionary["Unmap + Filt"] = splitWaltLine(line)[1]
    return value_dictionary
    
def extractHairpin(file_location):
    fd = open(file_location, 'r')
    value_dictionary = {}
    fields = ["Hairpin_Hits", "Test_Ambig", "Test_Correct", "Test_Incorrect", "Test_NotFound", "Test_UnmapFilt"] 
    for line in fd:
        for field in fields:
            if line.startswith(field):
                value_dictionary[field] = line.strip().split("\t")[1:3]
                break
    return value_dictionary
     
     
     
     
     
