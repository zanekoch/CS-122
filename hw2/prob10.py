import numpy as np
import argparse


def populate_matrix(matrixFile):
    '''
    reads in blosum matrix txt file and sets key: value pair as tuple(i, j) = score
    '''
    lines = matrixFile.readlines()
    matrixFile.close()
    dictaa = {}
    aminoacidstring = lines[0]
    aminoacidstring = aminoacidstring.split()

    i = 1
    while i <= (len(lines)-1):
        row = lines[i]
        row = row.split()

        j = 1
        for character in row[1:25]:
            dictaa[aminoacidstring[i-1],aminoacidstring[j-1]] = character
            j+=1
        i+=1
    #print(dictaa)
    return(dictaa)

def alignment(s1, s2, matrix):
    #create matrix that is len(s1)xlen(s1)
    #fill with 0s
    outputMatrix = np.zeros((len(s1)+1, len(s2)+1))
    backTrackMatrix = np.zeros((len(s1)+1, len(s2)+1))
    for i in range(len(s1)+1):
        outputMatrix[i][0]=5*i
    for j in range(len(s2)+1):
        outputMatrix[0][j]=5*j

    for j in range(1, len(s2)+1): #i:0->len(s1)-1, s1 is vertical
        for i in range(1, len(s1)+1): #j:0->len(s1)-1, s2 is horizontal
            deletion = outputMatrix[i-1][j] + 5
            insertion = outputMatrix[i][j-1] + 5
            identity = outputMatrix[i-1][j-1] if s2[j-1] == s1[i-1] else np.inf
            substitution = outputMatrix[i-1][j-1] - int(matrix[(s1[i-1],s2[j-1])]) if s2[j-1] != s1[i-1] else np.inf
            #print(s1[i-1],s2[j-1], int(matrix[(s1[i-1],s2[j-1])]))
            outputMatrix[i,j] = min(identity, deletion, insertion, substitution)
            if min(identity, deletion, insertion, substitution) == deletion:
                backTrackMatrix[i, j] == 'd'
            elif min(identity, deletion, insertion, substitution) == insertion:
                backTrackMatrix[i, j] == 'in'
            elif min(identity, deletion, insertion, substitution) == identity:
                backTrackMatrix[i, j] == 'id'
            else:
                backTrackMatrix[i, j] == 's'
    print(outputMatrix)
    return outputMatrix, backTrackMatrix

def back(backtrack, s1, s2, row, column):


    if row == 0 or column == 0:
        return
    if backtrack[row, column] == 'd':
        print('-')
        back(backtrack, s1, s2, row - 1, column)
    elif backtrack[row, column] == 'in':
        print(s1[row])
        back(backtrack, s1, s2, row, column - 1)
    elif backtrack[row, column] == 'id':
        print(s1[row])
        back(backtrack, s1, s2, row - 1, column - 1)
    elif backtrack[row, column] == 's':
        print(s1[row])
        back(backtrack, s1, s2, row - 1, column - 1)





def backtrace(s1, s2, matrix):
    current_row = len(s1) - 1
    current_column = len(s2) - 1
    changes = []
    while current_row > 0 or current_column > 0:
        if current_row == 0:
            pvs_row = -np.inf
        else:
            pvs_row = current_row - 1

        if current_column == 0:
            pvs_column = -np.inf
        else:
            pvs_column = current_column - 1

        try:
            insertion_dist = matrix[current_row, pvs_column]
        except IndexError:
            insertion_dist = np.inf

        try:
            deletion_dist = matrix[pvs_row, current_column]
        except IndexError:
            deletion_dist = np.inf

        try:
            if s1[current_row] == s2[current_column]:
                identity_dist = matrix[pvs_row, pvs_column]
            else:
                identity_dist = np.inf

            if s1[current_row] != s2[current_column]:
                substitution_dist = matrix[pvs_row, pvs_column]
            else:
                substitution_dist = np.inf
        except (TypeError, IndexError) as e:
            identity_dist = np.inf
            substitution_dist = np.inf
        min_dist = min(insertion_dist, deletion_dist, identity_dist, substitution_dist)

        s1_index = current_row
        if min_dist == identity_dist:
            print(s1[s1_index])
            current_row = pvs_row
            current_column = pvs_column
        elif min_dist == substitution_dist:
            changes.append(['SNP', s1[current_row], s2[current_column], s1_index])
            current_row = pvs_row
            current_column = pvs_column
        elif min_dist == insertion_dist:
            if len(changes) > 0 and changes[-1][0] == 'INS' and changes[-1][-1] == s1_index + 1:
                changes[-1][1] = s2[current_column] + changes[-1][1]
            else:
                changes.append(['INS', s2[current_column], s1_index + 1])
            current_column = pvs_column
        elif min_dist == deletion_dist:
            if len(changes) > 0 and changes[-1][0] == 'DEL' and changes[-1][-1] == s1_index + 1:
                changes[-1] = ['DEL', s1[current_row] + changes[-1][1], s1_index]
            else:
                changes.append(['DEL', s1[current_row], s1_index])
            current_row = pvs_row
        else:
            raise ValueError
        changes = sorted(changes, key=lambda change: change[-1])
        #print(str(changes))
        return changes



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', dest = 'matrixFile')
    args = parser.parse_args()
    matrixFile = args.matrixFile

    with open(args.matrixFile, 'r') as matrixFile:
        blosum_dict = populate_matrix(matrixFile)
    '''
    s1 = "ILYPRQSMICMSFCFWDMWKKDVPVVLMMFLERRQMQSVFSWLVTVKTDCGKGIYNHRKYLGLPTMTAGDWHWIKKQNDPHEWFQGRLETAWLHSTFLYWKYFECDAVKVCMDTFGLFGHCDWDQQIHTCTHENEPAIAFLDLYCRHSPMCDKLYPVWDMACQTCHFHHSWFCRNQEMWMKGDVDDWQWGYHYHTINSAQCNQWFKEICKDMGWDSVFPPRHNCQRHKKCMPALYAGIWMATDHACTFMVRLIYTENIAEWHQVYCYRSMNMFTCGNVCLRCKSWIFVKNYMMAPVVNDPMIEAFYKRCCILGKAWYDMWGICPVERKSHWEIYAKDLLSFESCCSQKKQNCYTDNWGLEYRLFFQSIQMNTDPHYCQTHVCWISAMFPIYSPFYTSGPKEFYMWLQARIDQNMHGHANHYVTSGNWDSVYTPEKRAGVFPVVVPVWYPPQMCNDYIKLTYECERFHVEGTFGCNRWDLGCRRYIIFQCPYCDTMKICYVDQWRSIKEGQFRMSGYPNHGYWFVHDDHTNEWCNQPVLAKFVRSKIVAICKKSQTVFHYAYTPGYNATWPQTNVCERMYGPHDNLLNNQQNVTFWWKMVPNCGMQILISCHNKMKWPTSHYVFMRLKCMHVLMQMEYLDHFTGPGEGDFCRNMQPYMHQDLHWEGSMRAILEYQAEHHRRAFRAELCAQYDQEIILWSGGWGVQDCGFHANYDGSLQVVSGEPCSMWCTTVMQYYADCWEKCMFA"
    s2 = "ILIPRQQMGCFPFPWHFDFCFWSAHHSLVVPLNPQMQTVFQNRGLDRVTVKTDCHDHRWKWIYNLGLPTMTAGDWHFIKKHVVRANNPHQWFQGRLTTAWLHSTFLYKKTEYCLVRHSNCCHCDWDQIIHTCAFIAFLDLYQRHWPMCDKLYCHFHHSWFCRNQEMSMDWNQWFPWDSVPRANCLEEGALIALYAGIWANSMKRDMKTDHACTVRLIYVCELHAWLKYCYTSINMLCGNVCLRCKSWIFVKLFYMYAPVVNTIEANSPHYYKRCCILGQGICPVERKSHCEIYAKDLLSFESCCSQKQNCYTDNWGLEYRLFFQHIQMECTDPHANRGWTSCQTAKYWHFNLDDRPPKEFYMWLQATPTDLCMYQHCLMFKIVKQNFRKQHGHANPAASTSGNWDSVYTPEKMAYKDWYVSHPPVDMRRNGSKMVPVWYPPGIWHWKQSYKLTYECFFTVPGRFHVEGTFGCNRWDHQPGTRRDRQANHQFQCPYSDTMAIWEHAYTYVDQWRSIKEGQMPMSGYPNHGQWNVHDDHTNEQERSPICNQPVLAKFVRSKNVSNHEICKKSQTVFHWACEAQTNVCERMLNNQHVAVKRNVTFWWQMVPNCLWSCHNKMTWPTRPEQHRLFFVKMRLKCMHEYLDVAPSDFCRNMQAYMHSMRAILEYQADFDLKRRLRAIAPMDLCAQYDQEIILWSGGYIYDQSLQVVSCEGCSYYADCYVKCINVKEKCMFA"
    '''
    s1 = "PLEASANTLY"
    s2 = "MEANLY"
    row = len(s1) - 1
    column = len(s2) - 1
    matrix, backtrack = alignment(s1, s2, blosum_dict)
    print(backtrack)
    changes = back(backtrack, s1, s2, row, column)
    print(changes)
