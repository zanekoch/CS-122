from constructBWT import *
from collections import defaultdict
import argparse
import time

def findNthOccurence(string, letter, n):
    start = string.find(letter)
    while start >=0 and n>1:
        start = string.find(letter, start+1)
        n-=1
    return start

def lastFirst(lastCol, i, firstCol):
    #find the "number letter" it is
    #occurence = occurenceNum(lastCol, lastCol[i], i)
    occurence = lastCol.count(lastCol[i], 0, i+1)
    location = findNthOccurence(firstCol, lastCol[i], occurence)
    return location

def BWMatching(lastCol, pattern):
    top = 0
    bottom = len(lastCol)-1
    while top <= bottom:
        if pattern != '': #does not enter if statement if pattern is empty
            symbol = pattern[len(pattern)-1] #last letter in pattern
            pattern = pattern[:-1] #remove last letter of pattern

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


if __name__ == '__main__':
    parser =<



    ref = "CCCAAAACACCAAATGTAGTAA"
    reads = ["CCC", "AAA"]
    ref += '$'

    lastCol, countArr, firstOccurrence, numList = globalizer(ref)


    locations = defaultdict(str)
    for read in reads:
        locations[read] = (BWMatching(lastCol, read))

    print(locations)
