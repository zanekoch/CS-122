

'''
state: sub problem is the largest subsequence of these two strings but shorter
relation: longest subseq between s1 and s2 (length i) is the longest subseq between s1[0:i] and s2[0:i] plus longest subseq between s1[i] s2[i]

btwn s1=A and s2=A is 1,A
btwn s1=AA s2=AC is 1,A => longest is longest between prev and if new letter is same plus that, if new letter is not same check if it creates a new longest
btwn s1=AAC s2=ACA is 2,AC or 2,AA =>longest is logest
btwn s1=AACC s2=ACAC is 3,ACC or 3, AAC
btwn s1=AACCT s2=ACACT is 4,ACCT or 4,AACT
btwn s1=AACCTT s2=ACACTG is 4,ACCT or 4,AACT
btwn s1=AACCTTG s2=ACACTGT is 5,ACCTT or 5,ACCTG or 5,AACTT or 5,AACTG
btwn s1=AACCTTGG s2=ACACTGTG is 6,ACCTTG or 6,ACCTGG or 6,AACTTG or 6,AACTGG
btwen s1=AACCTTGG s2=ACACTGTGA is 6,ACCTTG or 6,ACCTGG or 6,AACTTG or 6,AACTGG

longest between s1[0...m-1] s2[0...n-1] is:
    if(s1[m-1] == s2[n-1]):
        return s1[m-1] + longest between s1[0...m-2] s2[0...n-2]
    else:
        return max(longest(s1[0...m-2], s2[0...n-1]), s1[0...m-1] s2[0...n-2])
'''
'''
recursive way to do it
def L(s1, s2, n, m):
    if (n==0 or m==0):
        return ''
    if s1[n-1] == s2[m-1]:
        return L(s1[:-1], s2[:-1], n-1, m-1) + s1[n-1]
    else:
        return max(L(s1[:-1], s2, n-1, m), L(s1, s2[:-1], n, m-1))
'''
# Dynamic programming implementation of LCS problem

# Returns length of LCS for X[0..m-1], Y[0..n-1]
def lcs(X, Y, m, n):
    #build an [0...m] x [0...n] array of 0s
    L = [[0 for x in range(n+1)] for x in range(m+1)]
    # Following steps build L[m+1][n+1] in bottom up fashion. Note
    # that L[i][j] contains length of LCS of X[0..i-1] and Y[0..j-1]
    for i in range(m+1):
        for j in range(n+1):
            if i == 0 or j == 0: #if either string has run out of letters they cannot match so 0
                L[i][j] = 0
            elif X[i-1] == Y[j-1]: #if the letters match then the new longest is old longest plus new letter, so  one longer
                L[i][j] = L[i-1][j-1] + 1
            else:
                L[i][j] = max(L[i-1][j], L[i][j-1])#if the letters do not match try again including each letter alone
    # Following code is used to print LCS
    index = L[m][n] #set
    # Create a character array to store the lcs string
    lcs = [""] * (index+1)
    lcs[index] = ""
    # Start from the right-most-bottom-most corner and
    # one by one store characters in lcs[]
    i = m
    j = n
    while i > 0 and j > 0:
        # If current character in X[] and Y are same, then
        # current character is part of LCS
        if X[i-1] == Y[j-1]:
            lcs[index-1] = X[i-1]
            i-=1
            j-=1
            index-=1
        # If not same, then find the larger of two and
        # go in the direction of larger value
        elif L[i-1][j] > L[i][j-1]:
            i-=1
        else:
            j-=1
    print (''.join(lcs))


if __name__ == '__main__':
    X = "AGGAAGTTGGAGCCAACTTGGCATGCTGCGTTAGACGGGCAATGGACAACTTTGGGTAGACTCGTTTCATCAATTGTAACTGCCACGCCACCTGCGGAAGGGAGAGGGTCATGACTTCGGAAGCCACTGTTCCCGCATATTAACATTATGACGCGCGCAGCGTGAAATCGCAGAAAAAAGCAGCGGTAAAAGGCATGCTATTGCGCCGATAAATGTACTCTTGTTGTCGGCGCCTTTGTTTGATGAACGGAAGATCCGTGGGTCTGAACCGCTGATAGATATAAACCCGTTGGACCTCAGGAACGTATTTATCCACCCGTTGAGGGTTTCCCTGCAAAATAGCTTTCCTTTGGACCCCTCTTTAGCTGCATCTAGATCAGTCCCTCTATTTTCCCCCTTAAACGGATGTTATTTGACTTGCTGATGCTAATCCCCTATCCCATAAGCTCATATGCGAGGGAGGTGGGTTACGAAACGGTACCGTTATTGCCCCATTTCTGAACCGTCATAGCCATCTACCGCACAAGCGAGGCGGTTCAGTGTGGAATTTGGGGTTGAGGCCACGTACCTTCTTGCCGACTCGTACCTCGGAAATATGAAGATACCCACTGGAGTAGATACCCGTGTCGTTGTCCAAGAAAACGTCGACCTGTCTGGAAAATTTAACTTCAACCCGGATATGTACCCGATAGAAAATTCTACAACGCGCCGTGACCCGTAACCGCGTCCATTTCATCGGCACAACTTTCCTACTGCTACATTAGTTCCCAGAGGCAGACTTCCTTATACGAGCAAGACGCATCTTT"
    Y = "CACAGGTAATCTGGGCACACGCTGATTTTCTTACCTACCCGCCCAAGGCTCGGTGTGTAACTAAGAGGTCAACCGGTAGTCGCAATCACGGGCAGGTGATGCCGCCCGGGCAACACATGGGATCAAGTTAAATCAGTCATATGTGACTAACCTGATTCCAAACCCACATTCGGTATCCTGTGCAAAGGCATCCCGAACGCCCTCGCGACGACGCCCCTCAATCAATGATAATTGCAAAGACCAACGGATATAGTGTGTTTTTCCCCGCCCCTGAAAAGAACAAGTCGTACTCCCCGCTGTGATCGAATACGCGACCCGAGGTGATGCCCGCCGATTAGATTACGAGAGACAAAAACTTGTTGCCGGGCAACATACTTAAGAAGTTAAAAAGCAGCCCGCTACGTCCGGGTGGCGCGTTTGATATCCCTAGCAGCTGATATCCGTTGAGTAGCTAGGTCTGAAGGCATACTCAGATGGAGTGAGCGACGAATGAAGTCGGCAAATCGTATTGGTAATGCTCGCTGTATTTGCGGAAACTCAGTCACCGTTCGTTCCTGCGCGAAGTAGAAAGACAATACCACAGGTTGATTGCGGCCCAGTTACACGGACGGCGGGTCTAGAGACGTTGGACGTTGGCTACATTATGAACTAGGTCCTGCATGTTGGCTCATGTATGACACGTATTGTGACGGCTTAATCTGGTTAATTACTATGTTAAGAGATCTTCGGGATGGCCTTAGTGTAGAGATTACCTCAGGCAGCGAAAGAGACGCCGACGACTAGGAACGGAGGACTGACTTTCAGCTAGAGAATCACCAACTAACCTGGGTAATCACCCGGAACTTGGCTTATGTGGAGTACGGAAGATTCAAGAAGATCGACCGACCTATGGTTCCTCATTGC"
    m = len(X)
    n = len(Y)
    lcs(X, Y, m, n)
