"""
TO FIX:
dictionary alignment alg does not fit every read to genome because some reads' substrings of 10bp (out of 50bp) are not found in preproc reference genome dictionary
    -all 5 substrings, forward and backward, are not found for some reason
    -this is true for 528 reads
    -~6000 reads, so about 10% not alligned to genome at all
reason: these reads have significant mutations in them?


"""

import sys
import os
import argparse
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))
import numpy as np
from os.path import join
import time
from CM122_starter_code.helpers import read_reads, read_reference, pretty_print_aligned_reads_with_ref
from collections import defaultdict

#given a string (donorGenome) returns a dictionary of lists where the keys are the strings of length 10 and the strings are filled with where they occur
def preProcessRef(donorGenome):
    uniqueStrings = defaultdict(list)
    i = 0
    for i in range(len(donorGenome)-9):
        currentString = ""
        for j in range(10):
            currentString+=donorGenome[i+j]
        uniqueStrings[currentString].append(i)
    print(uniqueStrings)
    return uniqueStrings

#given a string returns where in the donor genome it is
def findString(string, uniqueStringsList):
    locations = uniqueStringsList.get(string)
    return(locations)

#returns a list of a string split into 5 pieces
def splitString(string):
    parts = [string[i:i+10] for i in range(0, len(string), 10)]
    return parts


#given a whole read (50bp), whole reference genome string, and where the read starts in ref genome, find the number of differences between them
def difference(read, ref, start):
    diffCounter = 0
    for i in range(3):
        if (read[i] == ref[i+start]):
            continue
        else:
            diffCounter+=1
    return diffCounter

#detect indels
#for the location and orientation where the minimum difference was found:
#       check if this start location appears in first subReads list of perfect matches and check if the end location +- some bound appears in last subreads readLocations
#if the distance between start of first subreads and start of second subread is not 30bp then there is an insertion or deletions
def indel_finder(minDiffLocation, revMinDiffLocation, uniqueStringsList, subStrings, subStringsRev, reversed):
    indelLocations = []
    indel = -10
    if reversed == False:
        endReadLoc = minDiffLocation + 30
        endLocs = [(endReadLoc+i) for i in (-5, -4, -3 , -2 , -1 , 1, 2, 3, 4, 5)]
        stringLocs = findString(subStrings[4], uniqueStringsList)
        if stringLocs == None:
            return [-1,-1]
        for loc in stringLocs:
            for l in endLocs:
                if loc == l:
                    indel = l
                    break
        if indel != -10:
            info = [minDiffLocation, indel]
        else:
            info =[-1,-1]
    if reversed == True:
        endReadLoc = revMinDiffLocation + 30
        stringLocs = findString(subStringsRev[4], uniqueStringsList)
        if stringLocs == None:
            return [-1,-1]
        endLocs = [(endReadLoc+i) for i in (-5, -4, -3 , -2 , -1 , 1, 2, 3, 4, 5)]
        for loc in stringLocs:
            for l in endLocs:
                if loc == l:
                    indel = l
                    break
        if indel != -10:
            info = [revMinDiffLocation, indel]
        else:
            info =[-1,-1]
    return info

#given a read:
#1) break it into 10 length pieces labeled 0,1,2,3,
#2) search for where each piece is in the donor genome, both forward and backward
#3) at each place where a piece fits perfectly align the whole read and record number of differences from genome, forward or backward depending on how piece fit
#4) to tie break where each fits check if its paired_end_reads also fits where it should before/after
def align(paired_end_reads, uniqueStringsList, ref):
    """
    :param paired_end_reads is a list of 2 element lists of two 50 character string
    :param uniqueStringsList is the preProcessRef dictionary
    :param ref is the reference genome as a string
    :return 1) a list of 2 element lists of the positoin where each paired read fits best_reads
            2) a list of 2 element lists of the 2 paired reads in their correct orientation
    """
    readLocations = []
    orientedReads = []
    unfoundReads = []
    for read_pair in paired_end_reads:
        onePairOfReadLocations = []
        onePairOfOrientedReads = []
        for read in read_pair:
            reversed = False
            subStrings = splitString(read) #list of 5, 10bp subStrings
            subNumber = 0
            minDiff = 1000
            minDiffLocation = -1
            unfound = 0
            for sub in subStrings: #each 10bp string in 50bp read
                foundStrings = findString(sub, uniqueStringsList) #list of where sub appears in reference
                if foundStrings ==  None: #case where substring has a mutation in it
                    unfound +=1
                    subNumber +=1
                    continue
                for subLoc in foundStrings: #check if the other reads around it match to the reference where it was foundStrings
                    readLocation = (subLoc - (10 * subNumber))
                    diff = difference(read, ref, readLocation)
                    if(diff < minDiff):
                        minDiff = diff
                        minDiffLocation = readLocation
                subNumber+=1
            #if (unfound == 5):
                #print(read)
            #for checking reversed reads
            reversed_read = read[::-1]
            subStringsRev = splitString(reversed_read)
            subNumber = 0
            revMinDiff = 1000
            revMinDiffLocation = -1
            for sub in subStringsRev: #each 10bp string in 50bp read
                foundStringsRev = findString(sub, uniqueStringsList) #list of where sub appears in reference
                if foundStringsRev ==  None: #case where substring has a mutation in it
                    unfound +=1
                    subNumber +=1
                    continue
                for subLoc in foundStringsRev: #check if the other reads around it match to the reference where it was foundStrings
                    readLocation = (subLoc - (10 * subNumber))
                    diff = difference(read, ref, readLocation)
                    if(diff < revMinDiff):
                        revMinDiff = diff
                        revMinDiffLocation = readLocation
                subNumber+=1

            if minDiff <= revMinDiff :
                onePairOfReadLocations.append(minDiffLocation)
                onePairOfOrientedReads.append(read)
                reversed = True
            else:
                onePairOfReadLocations.append(revMinDiffLocation)
                onePairOfOrientedReads.append(reversed_read)
                reversed = False
            if (unfound == 10):
                unfoundReads.append(read)

            #2 element list [start of read w indel, length of indel + for in - for del] --if no indel then [-1,-1]
            indels = (indel_finder(minDiffLocation, revMinDiffLocation, uniqueStringsList, subStrings, subStringsRev, reversed))
            #print (indels)

        orientedReads.append(onePairOfOrientedReads)
        readLocations.append(onePairOfReadLocations)
    return readLocations, orientedReads




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

    input_reads = read_reads(reads_fn)
    # This will take a while; you can use an array slice for example:
    #
    #   input_reads = reads[:300]
    #
    # to generate some data quickly.
    reference = read_reference(reference_fn)
    preProcessedReference = preProcessRef(reference) ##
    alignments, reads = align(input_reads, preProcessedReference, reference) ##

    output_str = pretty_print_aligned_reads_with_ref(reads, alignments, reference)
    with(open(output_fn, 'w')) as output_file:
        output_file.write(output_str)
