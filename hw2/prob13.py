def findNthOccurence(string, letter, n):
    start = string.find(letter)
    while start >=0 and n>1:
        start = string.find(letter, start+1)
        n-=1
    return start

def occurenceNum(string, letter, index):
    i = 0
    occurenceNum = 0
    while i <= index:
        if string[i] == letter:
            occurenceNum +=1
        i+=1
    return occurenceNum

#given a BWT'ed string, return the orignal string
def revBWT(text):
    #text is last col in cycled and sorted matrix
    lastCol = text
    #first column is sorted input
    firstCol = ''.join(sorted(lastCol))

    stopper = 0
    output = ''
    occurence = -1
    while(stopper != len(text)):
        if stopper == 0:
            dolLoc = findNthOccurence(lastCol, '$', 1)
            output += firstCol[dolLoc]
            occurence = occurenceNum(firstCol, firstCol[dolLoc], dolLoc)
            stopper +=1
        else:
            loc = findNthOccurence(lastCol, output[stopper-1], occurence)
            output += firstCol[loc]
            occurence = occurenceNum(firstCol, firstCol[loc], loc)
            stopper +=1
    return output



if __name__ == '__main__':
    text = "ACCTAGGCCCGGAGGACTGTAGTTCGAAGCACCAAACTACCTCCGAATGCTCACCCGGTATTACGTGACGGTATAAAGTACGAG$CCGAAGCGTTATTAAGGACCCAAAGTCGGTTAACGGCAACAGCAATACCGAGCTAGCCCTTTCGAGGGTAACTTGAAGTCATCGGGAGTACGAAGACAATGGAATAGCAATGATAGCGGGGATCTTTCTTCATCTAACACGCGATTCGACTTTTCGGGGGGCCCGGGCATCTTAGTTCAGCCGTAACCAGACATTGTGCAGGCGGCTGGCCGCACGTGTGCCAATAACCGCTACAATCGTCGTTGGCTCTCGGATTAACCGTCGATTAGCTGATCCCACAAGCCCCAGGGACTCCCTTGGGACCAATGTAGTTGAAGGCACCCTCTAGAACACGTACGAAATATCATGGGCCCCGCGTGAGGCGTCATGTGGTTGAATGTCCACCAATTTTCAGTACTCAGTAATACTGGTCCTCTCTCCTGATAGACCATAATTGTCGGTTTCACGCCAGTCGCTCCTCCGAGCAATGTCGTCGATTTGTGGTATAAGTAAAGTAGCTACTGGCGGCTTTAGCAGCTAAACTATCGAAGAAGACCCTCATCATAATCCACAACTACGGTAGAATTTTCACCTGAGCTGTTTATACAAGGACATCAGTTCCATATGACGTGCCGCGGTTGATGCACAGGAGCTTTGCGCTCCTAAATTATTTAGGCTAACTCATTCCCATGGCGGTATTGGCGCAGGTTCTTCCGACCAAAACTGGAACATCTGGGTAAGATCCCCTGGTTTCTACTGACCTTTCGCGATCGACGAATGCCCGTGTAGAGCTTACCCCAGACTATCACTGTAAAGTCTTTTGATCTTTATCCGCCGACGAATCCATATCCGCGTGAGAGATACACAGGAGTTGCGTATTAAGCTTCGCGACGTATAGTACATAGAATAGGACTAGTGCG"
    print(revBWT(text))
