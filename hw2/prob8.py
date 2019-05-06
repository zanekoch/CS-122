



def pathTo(n, m, down, right):
    """
    :down is a n x m+1 matrix with down[i,j] being the length of the downward street that can be walked from that position
    :right is a n+1 x m matrix with right[i,j] being the length that can be walked right from that point
    """
    dic = {(0,0):0}

    for i in range(1, n+1):
        dic[(i,0)]=dic[i-1,0]+down[i-1][0]
    for j in range(1,m+1):
        dic[(0,j)] = dic[(0,j-1)] + right[0][j-1]

    for i in range(1, n+1):
        for j in range(1, m+1):
            dic[(i,j)] = max(dic[(i-1,j)]+down[i-1][j], dic[(i,j-1)] + right[i][j-1])

    return dic[(n,m)]



if __name__ == "__main__":
    # read file and get parameters
    fin = open('rosalind_ba5b.txt','r')
    n = int(fin.readline().rstrip('\n') )
    m = int(fin.readline().rstrip('\n') )
    line = fin.readline().rstrip('\n')
    down_matix =[]
    while line != '-':
        down_matix.append([int(i) for i in line.split()])
        line = fin.readline().rstrip('\n')

    right_matrix = [[int(i) for i in line.split()] for line in fin.readlines()]

    output = pathTo(n,m,down_matix,right_matrix)

    print(output)

    """
    else:
        if right[i][j] + pathTo(i,j-1,n,m) >= down[i][j] + pathTo(i-1,j ,n,m):
             return = right[i][j] + max(pathTo(i,j+1,n,m), pathTo(i+1,j+1 ,n,m))
        else:
            return = down[i][j] + max(pathTo(i,j+1 ,n,m), pathTo(i+1,j+1,n,m))
    """
