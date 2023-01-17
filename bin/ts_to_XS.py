import sys
import itertools
import re

with open(sys.argv[1], "r") as myfile:
    reads = myfile.readlines()
    for read in reads:
        read = read.rstrip()
        tag = read.split(sep="\t")[1]
        if tag == "16":
            if "ts:A:+" in read:
                read = read.replace("ts:A:+", "XS:A:-")
            elif "ts:A:-" in read:
                read = read.replace("ts:A:-", "XS:A:+")
        elif tag == "0":
            read = read.replace("ts:A:", "XS:A:")
        print(read)
