from constructBWT import *
from collections import defaultdict
import argparse
import time
import sys
import os
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))
from helpers import read_reads, read_reference, pretty_print_aligned_reads_with_ref

#given a read fragment return a list of locations in matches perfectly in the genome
def BWMatching(readFrag):
    top = 0
    bottom = len(lastCol)-1
    while top <= bottom:
        if readFrag != '': #does not enter if statement if readFrag is empty
            symbol = readFrag[len(readFrag)-1] #last letter in readFrag
            readFrag = readFrag[:-1] #remove last letter of readFrag

            countArrCol = - 1
            if symbol == 'A':
                countArrCol = 1
            elif symbol == 'T':
                countArrCol = 2
            elif symbol == 'G':
                countArrCol = 3
            elif symbol == 'C':
                countArrCol = 4
            else: #$
                countArrCol = 0

            top = firstOccurrence[symbol] + countArr[top][countArrCol]
            bottom = firstOccurrence[symbol] + countArr[bottom + 1][countArrCol] - 1
        else:
            locations = [numList[i] for i in range(top, bottom+1)]
            #return bottom - top + 1,
            return locations

def splitRead(read):
    fragments = [read[i:i+10] for i in (0, 9, 19,29,39)]
    return fragments

def readMap(read):
    startLocations = [0, 9, 19, 29, 39]
    minMis = 1000
    minMisLoc = - 1
    #do stuff to first read forwards
    readFrags = splitRead(read)
    #find where each subread maps to
    for i in range(len(readFrags)):
        locations = BWMatching(readFrags[i])
        #for each of these locations check how well the whole read maps there
        if locations == None:
            continue
        for loc in locations:
            refLocation = loc - startLocations[i]
            if refLocation + 50 <= len(reference):
                mismatches = 0
                for j in range(len(read)):
                    if read[j] != reference[j + refLocation]:
                        mismatches += 1
                if mismatches < minMis:
                    minMis = mismatches
                    minMisLoc = refLocation
    return minMis, minMisLoc


def align(pairedReads):
    alignments = []
    oriented = []
    distancesBetween = []
    count = 0
    for readPair in pairedReads:
        onePairAlignments = []
        onePairOriented = []
        #find location each read in paired end reads fits best
        minMis, minMisLoc = readMap(readPair[0])
        revRead = readPair[1][::-1]
        minMisRev, minMisLocRev = readMap(revRead)

        '''if minMis > 20 or minMisRev > 20:
            continue'''
        #record places where reads do not match perfectly
        notPerfLocs = []
        if minMis > 0:
            notPerfLocs.append(minMisLoc)
        if minMisRev > 0:
            notPerfLocs.append(minMisLocRev)

        #record distance between for finding indels
        distanceBetweenReads = minMisLocRev - minMisLoc
        distancesBetween.append(distanceBetweenReads)

        #add to list of couples for input to pretty_print_aligned_reads_with_ref
        onePairAlignments.append(minMisLoc)
        onePairOriented.append(readPair[0])
        onePairAlignments.append(minMisLocRev)
        onePairOriented.append(revRead)
        alignments.append(onePairAlignments)
        oriented.append(onePairOriented)

        #timer stuff
        count += 1
        if count % 10000 == 0:
            time_passed = (time.clock())/60 #process_times
            print('{} reads aligned'.format(count), 'in {:.3} minutes'.format(time_passed))
            remaining_time = time_passed/count*(len(pairedReads)-count)
            print('Approximately {:.3} minutes remaining'.format(remaining_time))
    return alignments, oriented, distancesBetween, notPerfLocs


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
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

    #returns a list of couples of pairedendreads
    input_reads = read_reads(reads_fn)

    #returns a string
    reference = read_reference(reference_fn)
    reference += '$'

    #create BWT(reference) and auxillary data structures
    lastCol, countArr, firstOccurrence, numList = constructBWT(reference)

    #align the reads
    alignments, oriented, distancesBetween, notPerfLocs = align(input_reads)

    #do the pretty thing for pileup
    output_str = pretty_print_aligned_reads_with_ref(input_reads, alignments, reference)
    with(open(output_fn, 'w')) as output_file:
        output_file.write(output_str)

    '''
    reference = "CCCAAAACACCAAATGTAGTAA"
    reads = ["CCC", "AAA"]
    reference += '$'
    #create BWT(reference) and auxillary data structures
    lastCol, countArr, firstOccurrence, numList = globalizer(reference)
    locations = defaultdict(str)
    for read in reads:
        locations[read] = (BWMatching(lastCol, read))
    print(locations)
    '''
