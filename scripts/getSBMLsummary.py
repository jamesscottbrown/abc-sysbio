import abc
import sys


def main():
    
    source=[]
    for i in range(1,len(sys.argv)):
        source.append(sys.argv[i])

    abc.generateTemplate(source)

main()
