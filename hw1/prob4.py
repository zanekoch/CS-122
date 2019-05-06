import sys

def stripper(file_in, file_out):
    f = open(str(file_in), "r")
    fo = open(str(file_out), "w")
    for index, line in enumerate(f):
        if (index % 2 != 0):
            fo.write(line)
    return fo

stripper(sys.argv[1], sys.argv[2])
