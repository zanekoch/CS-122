import sys

def squarer(a, b):
    first = a*a
    second = b*b
    answer = first + second
    print(answer)
    return answer

squarer(int(sys.argv[1]), int(sys.argv[2]))
