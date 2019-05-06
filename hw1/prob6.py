import sys

def counter(s):
    dict={}
    for letter in s:
        if letter in dict.keys():
            dict[letter]+=1
        else:
            dict[letter]=1
    for letter, num in dict.items():
        print(letter, num)

counter(sys.argv[1])
