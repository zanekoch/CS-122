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

def BWMatching(lastCol, pattern, firstCol):
    top = 0
    bottom = len(lastCol)-1
    while top <= bottom:
        if pattern != '': #does not enter if statement if pattern is empty
            symbol = pattern[len(pattern)-1] #last letter in pattern
            pattern = pattern[:-1] #remove last letter of pattern
            #check if last column contains an instance of symbol
            if symbol in lastCol[top:bottom+1]:
                topCount = lastCol[0:top].count(symbol) + 1
                bottomCount = lastCol[0:bottom+1].count(symbol)
                topIndex = findNthOccurence(lastCol, symbol, topCount)
                bottomIndex = findNthOccurence(lastCol, symbol, bottomCount)
                top = lastFirst(lastCol, topIndex, firstCol)
                bottom = lastFirst(lastCol, bottomIndex, firstCol)
            else:
                return 0
        else:
            return bottom - top + 1


if __name__ == '__main__':
    #lastCol = "TCCTCTATGAGATCCTATTCTATGAAACCTTCA$GACCAAAATTCTCCGGC"
    lastCol = "ATTCTAATGATAGATAGTGCGTTGCACGTTTCACACATCCATAGTTTAGTCAGGCGCTCCCCATGCCGCGCGCTTTCTAATGCTCTTACCGGCTGTCGTTGTACCGTAGTCCGCGACCGGAGCTATGACTAGAGACTCTAACGCTTCCGACAAGGACGCCTCGGACTTTCTTCTGTGTTACGCTAATCGGTGTCTAGGACAGCCTCTGCGGAGGTCTATTAGGTGAATCGCTCGGTTCACGCGTATTCGAACAGTACGTGGATGCGCTGATTTTCGAGCAGGAGAGTTACTCTGTCGCACCGACTCCGCCGTATGCGAAGAGCCAGTCACACGGTCCAAGCGGAATGGATTGATGAACAACGGCGGGTCTGTAACACGTCGTTTCCGAGAAACAATTAACTGAGAAGGTACAAGCTTATTGTGACTGGAACCTACATATATTGCTCGTTGATCGGGGCCCCTAGCGTACAGAGGGGGGTCGACTGCCGACATCGCTTCCATATGTAGTACGGGAGCGACAGGCATTTTGCTCTTCGCCTCGCAACCAGCTCGATATGCTCAGCGTTCCCGCCGGTTAGAAATGAACTAGAGACCCTAGGCCTTAGCGCCCATGAATGCACTCCTCCACGCCAAGGTCTGAGGTGCTGACCGAGATCTGCGCGTTCTTGGGGTCAAAGACGGTCTAGTAGAGGATGAGGCGATTCATTTCTTAGCAATTGTCCCGTCCGCCAGGAGAAGTCCTTACACAACGTCTTTCACCATCATCTCACGTGAGCGAGTGTAGAGGACGCGAGTCTCAGCACGATGATCATCCGGGGCGATACGGCGCTCGGCGTTTTGCCAAGGCCAGGTGAGCTATTGAAGCAAAAAATGGGCATAATTCGCTGCATGTAATGGTTCGCTAAGCACAACACGTGGACGCTCGCAATACTGGGCCCAATCGTCGTTGGTGGCCCAGGAGCGGCGTTCCTCTCGTGAGGCGAAGAAGTGGTGCCGCCCGCGCGTCAGGATCGTAAAGACACGACCACTGATAAATGAAGCATTGGCGTGGCATACTTTTCTGATCATGGCCCTGGATGGATAGAGATCGAGAAATACGTATCCCGCCGCAGATGGACAGGCTAAACAGCTAGACGTAGTACGCAACAGGTTGTAACCGCTGCTTGTTGGCCGCAGTCGTAACTCGTAGATCATCAAGACGTTAACCGGGTCTTGTCAGGCTGGCAAGTCGTCTGGCGGTAGTGATCTCATCGAGGTTTGCGCTTCCAACGCTCCAGGTTCTACTAGACGGACACTTCCTCAACGCAGTATTTTTTCGCGTGTCTGGATGCCCCTCACTCGAAAGGGTTTGCAGATCAGTCGCGTCCAGACAACCTGCCAATGGATGAAAGATTAAGCCTTGATAACCAGACGGAGCACCACTCCTAACGAGGTGAGATTTACTCGAGCACAGATGAGACCTGCCGGAGCCATAACCGTCCCCTGTCATCGAGAGGATCCGGATCCGGCGGGTTACGATGCGTAAGTTAGTTAGATCCAAGACTTGTGATCAATCGAAATCGTAATTCCATACGGCCCCACCCGCGCGGCGGAACAGCCGAGAGGAATGCTAATGATTTGTGTATCATAGCGTGGTCACCGCCGCGAGGAAGATTAGTTAGAAGCAGATCGCTCAGTAGAAAACACACACGCAGAAGGCGGTTCGGTGACTGTCCTAGGATGCCTTCCTTTCCATTGGTCCAGT$AACGGTGCAAAGCAATCTCGCGTAAAGTTTGAAGGAGTCTGTGCTGCCGACATAAGAAAGTGCAACAGGTCCCCCGCCATAAGGTTAATATCGTTAGTAGTTGCGGCCAACTGAGAATACTAGTCATCCAGCCATCACCGGAAAGCGGAATCTCCGATGAACGTCTCCGACAAACAGACCAGACGTGCCCTTACGCGGGACAATTGCGAGTTTTAATGCGTCCTTACGAATTGGTAGGAGTGCGCCGTGGACAACCACACTACTTGAATGGCTTGTCAAAGTATCCCCGTCATCAGTGTGAGTACCGTACTTGTTGCTACAGCGCGGGCGGTTTGGCGTGACCGTTATCGACGATAAGTAAGCTGGGGAGTCCCTTCACTGAGAGTTATATGTACTTAAACCTGGACTAGAGGAGGTTAGTTACTCTAGCTCGGGGTAAGCCCATATCATGAGGTGAGAATTGGCATTACAGACTACCCCGTACACATGATTTCGCTCTACCATTCGGAATAGCAGCCGGCAGGATGCCATCTTACGCCCGTGCAAAGAAGAAGGATAGTAAGACGAGTTACACTGCCCGGGCTCGGGACGTTGCATTGGGCTTGGGTTATCTTACGCGTTTGTCGATTTTCTCCCTGCGGAATACGAATTCGTTCTTCTGGACAATGTCTGTGGCAAATATGTTCTGACGCAGTTAAATTCTGCGAGGCCATAATACGGCCCAAGCTTATCCAACGCAGCAGCGTAGTGTAAGTCTGATCGAAGTGGTCATGAGACCTTCACGCGATCCTATTCAGGGTCTAAATGAATACCTACTAGACCGTCTGATATTCTGCTCCAGGGGGTCTCCTAGATATATTCGTGCCCTTGTATCCGCAACCTCCGAGGCTGGGCATGAGGACGACCGACCCTTGGATAACAGAGTGTGTCGCCGGTACGGAGAAGTTCGCTTAACGGCACTGCCACCTGGGTATATTAACAATGAGCTAGCACCCGCAGGTCATGTTTTACTCTGAGCGCGTGGACTTTGAAATTTAGGCTGCGACGTCTAGCCCGGGGGTTGCCCCGACTGCCCAGATCCTTTGAGCAATGGCCAGAGTGCACAGCGGAATGCCATTCGCCGCACTCTGCAGCTTGAGATAGTTCAGACTCCACGCGTTGCGGAATTAAGCGCGGATGCACCGCTCCACATTTGTAACTGTCTGTTGTAGAGTCGACGTAATATCTACAGAACAGGGGTACGATCGCGGAGGCGGGCGACTACGTGATCAGACAGGTGTTGGGCAACCGAACGTTCTGCATTTTATCCTGGTTGAGACAAACGAGACTCCGACATTCAGCAATCATGCTTGGGATCCAGCTTAATAGATGGTGCAAGTACAGCCCAAGCATAGCAGTCGGACGGGTCGGTGACGTACGAAGCTGCCTCAGGTAGTAACGCGTGCCGTTCACCACAGCAATTTCGTCTCAGCTCCCTGGACGCCCAGAAAGCCTCCTGTGATTAGCAGGGTTGGGCGACGTTCGAGAGTTTTACAGATATCTTCATGATTCTAGTGTTAAACAAGGCGGCCAAGCCAAGGGTTCCGGGCCCGCGGTACGCTCTCATTCTGTGTCGTCCACTTGCTATAGAAGTGAATACGCATGAGGCTTTTTGCATGTGGGGCCTACCGAATGGAGCTCCCATAGTTATAGTGCTACTAGACGATACTAGCCAAAAAGGACTACCGCTACTATCGTGCTTGGGGACGACAAAAACAATATCCTGCTTCCGTATCGCTGGGGGCGATCCTAGACCGACAGGCAGTTTAACTTTGCGACTGGGGAGCCATCCAATAGCCTACTTGTTACTCGGATTCAGTTGAACTATACACAGTATAGACCCATGGGCGATGGATCGCCCTGAACGACGGAGATACCCCTCGCTGGGTCTGGAATGGCTTAAATTACCCCACTACATCCAGCTCGGAGACTGGCAACTTATCAAACCGCGTGCAGCTCCGGTACCGAACTCCCGGGTCTGGGACGCGATCAAGTATGGATACTGGCTAGACGAGCCGAAAACTTAACGACGATAATTCCGTGTGCTTCATTCATATGCGTAGAAGCCGCATCTTCCGCATATTGGCGTAGTTACAACCAGCTTTAGTGTGGTCATGAAATCTTGCGAGGAACCTACGCAGTACGCTTAGGAAAAAGTTTCCGGCAGCGTTAAGATTGTACTGCAGGTCGAGCTAGCTCTCCCAACGCTTCTGCCACTGCGCACCCCGCGAAACTTACCTCGGCCGCGTTGCACTGTTTTGAATCACCGGTCGTTTAAGTGCCTGCGGCCATGTGGCTGAGGCATATCTGAGCTAAGGCTGTTGAGCCCAGAAGGCCTTTGACACAGGCCTGTCATATGAAAGATCCGTAACATAGGCCACGGGTATGCTCTTTTAACTGTAATTATGATGATTATATACACTTCCGAATACCTTATATGCTAAATCAAGAAGAGACGGGAAGCCTAGATATTAAATTTTGCAGACCTTCCTCCGGATCTTCAATCCGCCTAAAAGAGAATACGGGCCCTCCGTGTTCGGGTGTCATAAACTCTACCTTAAACGACTTTCTACTGTATCAAGTTGTCCAGACCGGTCCGAGATGATTCGTGACGTTCCCGAGACCACAGGAGTTAGCGCCAGGTAGTCAAAACTGGATCGTTCATGTAAGTGCGGTAGACCACACGCTTTATCCGGTTAGACTTTGACCGGACGCACTTCGTGGTCCGTGGCGTGCGTCTTAGGGCCAAACTAGGCCTAGTACAGTTACATTAGACATCACCAATAAGTGATTACCGTGTCACTACTATTTTAGATGCACATCGGGGTCAGCGTTTAGGATGTGGTATTCGGCACATACTTATAGGGTAACGCAATCGGTCGTGCATCCTAAGCGTCTGTGGTCAGAAAGCCATTGACGTAGTCGTTCTTAATTTCTCATACCTACAGGACACTAGAGAGCTACGGACCCCCCGAACCATCTCATAAGGGATATCTAAACTACCCTCGTGCACCGATGGGCTGTCCGATGAGCAAAAGAGTCCTCCAATATGGAGTTCTCCATGCGAGGTCTGGGGGAGACCGTAATAGATACTGCCCACCTGGCCTGTTATGGGATGAGAAATGATAATAGTGCCGTCGTCCGAACGCTCCAGGTACCGTAGAATTCAGGGTCCATAAGTCGCCTTAATCGAAGAAAGCCTCAAGTGAGCGTTAACTAGTGACGTCAGCGGTTCGCGGTGTTCTTGCTCACGGGGGATCGTCTTTGCGAGTCTATCCCGCGTCCTCATGTAATGTGACCGGAATGTGCAAGTTTTCTAATATACACAAGCTGTCTGAAATTTCTGTCTTACTGTTTACCCAATTTAACATTCTTAAAGCTTTGTTCTATATACAGACCGGGGTTGGCGTACGGTATCATTGTTAGGAACTAATTATGTAAACTGAGCTCCAAGACGGGCGGTTCGCCCGTGAAAGTCCAATGACGTCCGAGGAGGAGGCGCCAGTCTTCACACGGTTACTAGTTCTTTCTCCTCGTATTGCTTGCTTAGACGACATACCTTGCACTAGCCGGTGGGTCCCGCCTCCACCGAAATTTTGTTGCGCCGCATCTCCTAGTGCAGCAGGTCTACGCCTATCCGAGATGTTCGACAATATTCTGGTCTTAACACGTCTACGATCCAATCCGAAAATGGTGCTCCGCTGCACTGGCCTGCTATGGACGTGCGGAATGCTGCGCAGCCACAAAGGGTATGCCCACTGCGAACATCCGAGAACCCTTGGATCAACGTCACTCGGTCTGAGGACCACAGAGTAGTGTGATCTGGAGGCTTATCGCCTAAACCGGAATCCTGTGGGTTGTTCGCAAACGAAGGGCATGCGGAGTCCACAGTCCCGCGAACCCGCGTTATATGCAAAGTCGGAGCTAGCGATATCGGCCAAATCGGACAGTGTACCTCCCCATGGACCAGTGCGATAAGTCCAGTCCGTGCTATCTAACCGTAATTGGTTTGGGCCCTCTGTGATGCGTTTAAACAACTCGAGCAGTAGCAAATGACGATACGTAAAACATTGTCGATGAGTGACCTGGTGTCATCCGACGCTGTCTAAGGACTACACCGCGCCCCTATAGGAGATATATGTATTTCACTTCTCCTGGTCTCCATTGCTTAAGCAGTGGCCATGGTCGGCTACGTTGAGAATATTCTCGATCGATCG"
    firstCol = ''.join(sorted(lastCol))
    patterns = []
    with open("reads2.txt", "r") as readFile:
        patterns = readFile.read().split()
    #print(patterns)
    #patterns = ["CCT", "CAC","GAG","CAG","ATC"]
    data =[]
    for bruh in patterns:
        data.append(str(BWMatching(lastCol, bruh, firstCol)))
    print(*data, sep=' ')
