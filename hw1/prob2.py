import sys
def slicer(string, a, b, c, d):
    first = ""
    second = ""
    for letter in string[a:b+1:1]:
        first+=letter
    for letter in string[c:d+1:1]:
        second+=letter
    print(first," ",second)
slicer(sys.argv[1],int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]))
