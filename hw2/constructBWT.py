import numpy as np
from collections import defaultdict

#i is the index of the row
def keyCreator(i):
    #gets the cycled row
    row = text[i:] + text[:i]
    return row


def allInOne(i, countArr, numList, firstOccurrence):
    row = text[numList[i]:] + text[:numList[i]]

    #update countArr
    if row[len(text) - 1] == 'A':
        countArr[i+1][1] =  countArr[i][1] + 1
        countArr[i+1][0] =  countArr[i][0]
        countArr[i+1][2] =  countArr[i][2]
        countArr[i+1][3] =  countArr[i][3]
        countArr[i+1][4] =  countArr[i][4]
    elif row[len(text) - 1] == 'T':
        countArr[i+1][2] =  countArr[i][2] + 1
        countArr[i+1][0] =  countArr[i][0]
        countArr[i+1][1] =  countArr[i][1]
        countArr[i+1][3] =  countArr[i][3]
        countArr[i+1][4] =  countArr[i][4]
    elif row[len(text) - 1] == 'G':
        countArr[i+1][3] =  countArr[i][3] + 1
        countArr[i+1][0] =  countArr[i][0]
        countArr[i+1][1] =  countArr[i][1]
        countArr[i+1][2] =  countArr[i][2]
        countArr[i+1][4] =  countArr[i][4]
    elif row[len(text) - 1] == 'C':
        countArr[i+1][4] =  countArr[i][4] + 1
        countArr[i+1][0] =  countArr[i][0]
        countArr[i+1][1] =  countArr[i][1]
        countArr[i+1][3] =  countArr[i][3]
        countArr[i+1][2] =  countArr[i][2]
    else:
        countArr[i+1][0] =  countArr[i][0] + 1
        countArr[i+1][4] =  countArr[i][4]
        countArr[i+1][1] =  countArr[i][1]
        countArr[i+1][3] =  countArr[i][3]
        countArr[i+1][2] =  countArr[i][2]

    #update firstOccurrence
    if row[0] == 'A' and 'A' not in firstOccurrence:
        firstOccurrence['A'] = i
    if row[0] == 'C' and 'C' not in firstOccurrence:
        firstOccurrence['C'] = i
    if row[0] == 'T' and 'T' not in firstOccurrence:
        firstOccurrence['T'] = i
    if row[0] == 'G' and 'G' not in firstOccurrence:
        firstOccurrence['G'] = i
    if row[0] == '$' and '$' not in firstOccurrence:
        firstOccurrence['$'] = i


    #add letter in lastcol to lastCol
    return row[len(text) - 1]



def BWT(text):
    #list of numbers 0->len(text)-1 that correspond to the unsorted rows of cycled matrix
    #intiially index i corresponds to row number i
    numList = [i for i in range(0,len(text))]

    #sorts the list in place based on the cycled strings, now index i corresponds to some new row

    numList.sort(key=keyCreator)

    lastCol = ""

    #getting the last column of the matrix, creating count data structure, and first occurence data structure
    #count is a len(text)+1 x 5 matrix (one row taller than text)
    #firstOccurrence is a dictionary where the key is the letter and the value is the location where it first appears in the first col of BWT
    countArr = np.zeros(shape=(len(text)+1, 5), dtype=int)
    firstOccurrence = defaultdict(str)
    for i in range(len(numList)):
        lastCol += allInOne(i, countArr, numList, firstOccurrence)

    return lastCol, countArr, firstOccurrence, numList

def globalizer(textt):
    global text
    text = textt
    lastCol, countArr, firstOccurrence, numList = BWT(text)
    return lastCol, countArr, firstOccurrence, numList

if __name__ == '__main__':
    text = "GCGTGCCTGGTCA$"
    lastCol, countArr, firstOccurrence, numList = BWT(text)
    print(lastCol)
    print(countArr)
    print(firstOccurrence)
