import sys

def permuter(s, prefix, n, k):
    f = open("output.txt", "a")
    if k == 0:
        f.write(prefix + '\n')
        return
    for letter in s:
        newPrefix = ""
        newPrefix = prefix + letter
        permuter(s, newPrefix, n, k - 1)
empty = ""
permuter(sys.argv[1], empty, int(sys.argv[2]),int(sys.argv[2]))
