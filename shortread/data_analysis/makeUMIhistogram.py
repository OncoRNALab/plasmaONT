#!/usr/bin/env python
import pandas as pd
import os
from io import StringIO
import argparse
from gzip import open as gzopen
import time
from collections import Counter
start = time.time()
#inputfile = "/Users/rmvpaeme/test/UMI_trim/example/UMI_test.bam"
#outputfile = "/Users/rmvpaeme/test/UMI_trim/example/UMI_test_umtag.bam"
#UMIfile = "/Users/rmvpaeme/test/UMI_trim/UMI_list_RNA_exome.txt"
#inputfile = "/Users/rmvpaeme/Repos/ExtractUMI/example/UMI_test.fastq.gz"

## Read arguments from command line
parser = argparse.ArgumentParser(
"""Extracts UMIs from RNA exome data.

Makes a histogram from the first n basepairs in a fastq file
""")

parser.add_argument('-i', '--input', help = '.fastq.gz', required = True)
parser.add_argument('-o', '--output', help = '.txt output file', required = True)
parser.add_argument('-l', '--length', help = 'UMI length', required = True)
args = parser.parse_args()
inputfile = args.input
outputfile = args.output
length = int(args.length)

if not os.path.exists(inputfile):
    raise SystemError("Error: File does not exist\n")

def process(lines=None):
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}


print(
"""
makeUMIhistogram

Running in FASTQ mode.
Input: %s
Output: %s
""" % (inputfile, outputfile))

#out = gzip.open(outputfile, 'wt')
n = 4
lines = []
reads = 0
UMI_counts = []
for line in StringIO(gzopen(inputfile).read().decode("utf-8")):
    lines.append(line.rstrip())
    if len(lines) == n:
        record = process(lines)
        UMI = record["sequence"][0:length] 
        UMI_counts.append(UMI)
        lines = []
        reads = reads + 1
        if reads % 100000 == 0:
            print("%s reads processed ..." % (reads))


UMI_counts_df = pd.DataFrame.from_dict(dict(Counter(UMI_counts)), columns = ["counts"], orient = "index").sort_values(by = "counts", ascending = False)
UMI_counts_df.index.name = "UMI"
UMI_counts_df.to_csv(outputfile)


end = time.time()
time = end-start
print("Finished running. Runtime %.2f seconds." % time)
