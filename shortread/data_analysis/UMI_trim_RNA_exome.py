#!/usr/bin/env python
import pandas as pd
import numpy as np
import sys
import os
from Bio import SeqIO, bgzf
from io import StringIO
import argparse
import gzip
from gzip import open as gzopen
import pysam
import time
start = time.time()
#inputfile = "/Users/rmvpaeme/test/UMI_trim/example/UMI_test.bam"
#outputfile = "/Users/rmvpaeme/test/UMI_trim/example/UMI_test_umtag.bam"
#UMIfile = "/Users/rmvpaeme/test/UMI_trim/UMI_list_RNA_exome.txt"

## Read arguments from command line
parser = argparse.ArgumentParser(
"""Extracts UMIs from RNA exome data.

This program takes a list of UMIs as input and checks whether a read begins with any of those UMIs (FASTQ) or if the end of the read name contains a UMI sequence (BAM).

FASTQ-mode: The UMI is extracted from the read and added to the read name, and 8 nt are trimmed off.
BAM-mode: The UMI is extracted from the read name and added as a "UM" tag.

If a UMI is 6 nt, there is an N added to the UMI (so eventually all the UMIs are 7 nt long). Reads that do not begin with a UMI are discarded.
""")

parser.add_argument('-i', '--input', help = '.fastq.gz or sorted and indexed bam file', required = True)
parser.add_argument('-o', '--output', help = '.fastq.gz file (will be gzipped) or bam file', required = True)
parser.add_argument('-u', '--umifile', help = ".txt file with UMI sequence on each line", required = True)
parser.add_argument('-b', '--bam', help = "if input is bam file", action="store_true")
args = parser.parse_args()
inputfile = args.input
outputfile = args.output
UMIfile = args.umifile
bamfile = args.bam

if not os.path.exists(inputfile):
    raise SystemError("Error: File does not exist\n")

def process(lines=None):
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}

UMI_list = pd.read_csv(UMIfile , sep = ",", header = None)

if bamfile == False:
    print(
    """
ExtractUMI

    Running in FASTQ mode.
    Input: %s
    Output: %s
    UMI file: %s
    """ % (inputfile, outputfile, UMIfile))
    out = gzip.open(outputfile, 'wt')
    n = 4
    lines = []
    for line in StringIO(gzopen(inputfile).read().decode("utf-8")):
        lines.append(line.rstrip())
        if len(lines) == n:
            record = process(lines)
            for UMI in UMI_list[0]:
                if (len(UMI) == 6) & (str(record["sequence"]).startswith(UMI)):
                    out.write(record["name"]+ '_'+ UMI + 'N' +'\n')
                    out.write(record["sequence"][8:]+ '\n')
                    out.write(record["optional"]+ '\n')
                    out.write(record["quality"][8:]+ '\n')
                    lines = []
                    break
                elif (len(UMI) == 7) & (str(record["sequence"]).startswith(UMI)):
                    out.write(record["name"]+ '_'+ UMI +'\n')
                    out.write(record["sequence"][8:]+ '\n')
                    out.write(record["optional"]+ '\n')
                    out.write(record["quality"][8:]+ '\n')
                    lines = []
                    break
                else:
                    lines = []
                    break

    out.close()

else:
    print(
    """
ExtractUMI

    Running in BAM mode.
    Input: %s
    Output: %s
    UMI file: %s
    """ % (inputfile, outputfile, UMIfile))
    bam = pysam.AlignmentFile(inputfile)
    output_umi = pysam.AlignmentFile(outputfile, "wb", template=bam)
    for read in bam.fetch():
        UMI_from_read = str(read).split("\t")[0].split("_")[1]
        for UMI in UMI_list[0]:
            if (len(UMI) == 6) & (str(UMI_from_read).startswith(UMI)):
                read.tags += [('UM', UMI + 'N')]
                output_umi.write(read)
                break
            if (len(UMI) == 7) & (str(UMI_from_read).startswith(UMI)):
                read.tags += [('UM', UMI)]
                output_umi.write(read)
                break
    bam.close()
    output_umi.close()

end = time.time()
time = end-start
print("Finished running. Runtime %.2f seconds." % time)
