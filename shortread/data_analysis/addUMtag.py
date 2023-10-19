#!/usr/bin/env python
import sys
import os
import argparse
import pysam
import glob
import time

start = time.time()

## Read arguments from command line
parser = argparse.ArgumentParser(
"""Extracts UMI from read name and adds as UM-tag.
""")

parser.add_argument('-i', '--input', help = 'input', required = True)
parser.add_argument('-o', '--output', help = 'output', required = True)
args = parser.parse_args()
inputfile = args.input
outputfile = args.output

if not os.path.exists(inputfile):
    raise SystemError("Error: File does not exist\n")

print(
"""
extractUMtag
    Input: %s
    Output: %s
""" % (inputfile, outputfile))

import pysam
totalReads = 0
bam = pysam.AlignmentFile(inputfile)
output_umi = pysam.AlignmentFile(outputfile, "wb", template=bam)
for read in bam.fetch():
    UMI = str(read).split("\t")[0].split("_")[1]
    read.tags += [('UM', UMI)]
    output_umi.write(read)
    totalReads += 1
bam.close()
output_umi.close()

print("Number of reads in BAM file: %s " % totalReads)
end = time.time()
time = end-start
print("Finished running. Runtime %.2f seconds." % time)
