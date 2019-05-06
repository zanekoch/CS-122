import sys

def summer(a, b):
    answer = 0
    if(a % 2 ==0):
        for num in range(a+1,b,2):
            print(num)
            answer+=num
        answer+=b
    else:
        for num in range(a,b,2):
            print(num)
            answer+=num
        answer+=b
    print(answer)

summer(int(sys.argv[1]), int(sys.argv[2]))
