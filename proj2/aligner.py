import sys
import os
import argparse
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))
import numpy as np
from os.path import join
import time
from helpers import read_reads, read_reference, pretty_print_aligned_reads_with_ref
from collections import defaultdict

def preProc(ref):
    """
    :param ref is a string version of reference genome
    :returns a dictionary with keys that are 10bp strings that leads to lists of where this string occurs
    """
    uniquePieces = defaultdict(list)
    for i in range(len(ref)-9):
        currentString = ""
        for j in range(10):
            currentString+=ref[i+j]
        uniquePieces[currentString].append(i) #appends the location the string occurs at
    return uniquePieces


def findRead(read, preProcRef):
    """
    :param read a input string (subread)
    :param preProcRef the preproc dictionary of reference genome
    :return list of where in genome is is (list can be 0, 1, or many elements long)
    """
    locations = preProcRef.get(read)
    return(locations)


def splitRead(read):
    """
    :param read a input stirng
    :returns a list of 10 char long pieces of strings
    """
    parts = [read[i:i+10] for i in range(0, len(read), 10)]
    return parts

def readMap(read, ref, preProcRef):
    """
    :param read is a 50bp read
    :param ref is the reference genome string
    :param preProcRef is the dictionary of reference genome
    :return the minimum number of mismatches the read had and the location this happened
    """
    numberSubReads = 3
    subReadLocs = [int(i*len(read)/numberSubReads) for i in range(numberSubReads)]
    subReads = [read[loc:loc + 10] for loc in subReadLocs] #10 bc that is how long want reads to be
                                                           #makes 3 reads of length 10 starting at positions 0, 16 and 33 in read
    minMis = 1000
    minMisLoc = -1
    for i in range(numberSubReads):
        if subReads[i] in preProcRef:
            for loc in preProcRef[subReads[i]]:
                startOfWholeRead = loc - subReadLocs[i]
                if startOfWholeRead + len(read) <= len(ref):
                    mismatches = 0
                    for j in range(len(read)):
                        if read[j] != ref[j + startOfWholeRead]:
                            mismatches += 1
                    if mismatches < minMis:
                        minMis = mismatches
                        minMisLoc = startOfWholeRead
    return minMis, minMisLoc


def align(pairedReads, ref, preProcRef):
    """
    :param pairedReads is a list of elements lists of two 50 char strings
    :param preProcRef is a dictionary with keys that are 10bp strings that leads to lists of where this string occurs
    :param ref is reference genome string
    :return 1) a list of 2 element lists of the positoin where each paired read fits best_reads
            2) a list of 2 element lists of the 2 paired reads in their correct orientation
    """
    alignments = [] #list of two element lists of locations of reads
    oriented = [] #list of two element lists of reads in their correct orientations
    count = 0
    for readPair in pairedReads:
        onePairAlignments = []
        onePairOriented = []
        for read in readPair:
            minMis, minMisLoc = readMap(read, ref, preProcRef)
            revRead = read[::-1]
            minMisRev, minMisLocRev = readMap(revRead, ref, preProcRef)

            if minMisRev < minMis:
                minMis = minMisRev
                minMisLoc = minMisLocRev
                read = revRead

            if minMis <= 2: #if less than tolerance
                onePairAlignments.append(minMisLoc)
                onePairOriented.append(read)
        #check if both reads are there
        if len(onePairAlignments) < 2:
            continue
        else:
            alignments.append(onePairAlignments)
            oriented.append(onePairOriented)
            count += 1
            if count % 10000 == 0:
                time_passed = (time.clock())/60 #process_times
                print('{} reads aligned'.format(count), 'in {:.3} minutes'.format(time_passed))
                remaining_time = time_passed/count*(len(pairedReads)-count)
                print('Approximately {:.3} minutes remaining'.format(remaining_time))
    return alignments, oriented




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_aligner.py takes in data for homework assignment 1 consisting '
                                     'of a genome and a set of reads and aligns the reads to the reference genome.')
    parser.add_argument('-g', '--referenceGenome', required=True, dest='reference_file',
                        help='File containing a reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    args = parser.parse_args()
    reference_fn = args.reference_file
    reads_fn = args.reads_file
    output_fn = args.output_file

    input_reads = read_reads(reads_fn) #returns a list of couples of paired_end_reads
    # This will take a while; you can use an array slice for example:
    #
    #   input_reads = reads[:300]
    #
    # to generate some data quickly.

    reference = read_reference(reference_fn) #returns a string of straight up reference dawg

    preProc = preProc(reference)

    alignments, reads = align(input_reads, reference, preProc)

    output_str = pretty_print_aligned_reads_with_ref(reads, alignments, reference)
    with(open(output_fn, 'w')) as output_file:
        output_file.write(output_str)
