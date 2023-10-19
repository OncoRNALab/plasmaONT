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
import glob
import time
import itertools
import multiprocessing
from multiprocessing import Process, Manager, Pool
from collections import Counter
import subprocess
cpuCount = (multiprocessing.cpu_count() - 2)

start = time.time()
# inputfile = "/Users/rmvpaeme/Repos/ExtractUMI/example/UMI_test.bam"
# inputfile = "/Users/rmvpaeme/Repos/ExtractUMI/example/small_picoV3.bam"
# outputfolder = "/Users/rmvpaeme/Repos/ExtractUMI/example/".rstrip("/")
# UMIfile = "/Users/rmvpaeme/Repos/ExtractUMI/UMI_list_RNA_exome.txt"
#os.chdir("/Users/rmvpaeme/Repos/ExtractUMI")

## Read arguments from command line
parser = argparse.ArgumentParser(
"""Extracts UMIs from RNA exome UMI data and picoV3 libraries (not compatible with other UMI types).

In fixed mode, this program takes a list of UMIs as input and checks whether a read in bam file contains a UMI in a user-provided list. If there is no match, that read is discarded.
In random mode, this program removes all UMIs that contain "N".

If a UMI is 7 nt, they are trimmed to the first 6 nt.

The UMI reference file only needs to contain the R1 UMI. If the data is paired end, all combinations are automatically generated from that for PE reads (e.g. 120 UMIs in SE sequencing = 120^2 UMIs in PE sequencing). This requires of course that the UMIs (and length) in R1 and R2 are identical.

The program requires that the UMI is encoded in the read name (e.g. NS500361:616:HJ3LFAFXY:3:21412:2959:9891_TACCTGGAACCT), which is automatically done if the fastq was preprocessed with umitools.
""")

parser.add_argument('-i', '--input', help = 'sorted and indexed bam file', required = True)
parser.add_argument('-o', '--outputfolder', help = 'bam file after denoising', required = True)
parser.add_argument('-u', '--umifile', help = ".txt file with UMI sequence on each line (required if -m fixed)", required = False)
parser.add_argument('-m', '--mode', nargs=1, choices=['fixed', 'random'], help = "are the UMIs predetermined (fixed) or random?", required = False)
parser.add_argument('-p', '--paired', help = "if paired end reads", required = True, action='store_true')
parser.add_argument('-c', '--cores', help = "number of cores for parallel processing. Default = cpuCount - 2.", default = (multiprocessing.cpu_count() - 2))
parser.add_argument('-s', '--chunksize', help = "number of reads in every chunk. Default 1e6.", default = 1000000)
args = parser.parse_args()
inputfile = args.input
outputfolder = args.outputfolder.rstrip("/")
mode = args.mode
mode = "".join(mode)
if mode == "fixed":
    UMIfile = args.umifile
elif mode == "random":
    UMIfile = "running in random mode"
else:
    raise SystemError("specify mode\n")
pairedEnd = args.paired
cpuCount = int(args.cores)
chunk_size = int(args.chunksize)
logKept = outputfolder + "/" + inputfile.split('/')[-1].split('.')[0] + "_noNoise_keptUMIs.txt"
logRejected = outputfolder + "/" + inputfile.split('/')[-1].split('.')[0] + "_noNoise_removedUMIs.txt"
tmp_pattern = outputfolder + "/" + inputfile.split('/')[-1].split('.')[0] + "_tmp%d.bam"
outfile_pattern = outputfolder + "/" + inputfile.split('/')[-1].split('.')[0] + "_noNoise.bam"

if not os.path.exists(inputfile):
    raise SystemError("Error: File does not exist\n")

if mode == "fixed":
    UMI_list = pd.read_csv(UMIfile , sep = ",", header = None)
    UMI_list_6nt = UMI_list.apply(lambda x: x.str.slice(0,6), axis = 1)

    if pairedEnd == True:
        peUMI_list = []
        for p in itertools.product(UMI_list_6nt[0], UMI_list_6nt[0]):
            umi = ''.join(p)
            peUMI_list.append(umi)
        UMI_list_6nt = pd.DataFrame(peUMI_list)

#UMI_list_6nt.to_csv("/Users/rmvpaeme/Repos/ExtractUMI/UMI_list_RNA_exom_6nt.txt", index = False)

print(
"""
denoiseUMI
    Input: %s
    Output folder: %s
    Output file: %s
    Mode: %s
    UMI reference file: %s
    Log files:
        Kept UMIs: %s
        Removed UMIs: %s
    Paired-end: %s
    Cores: %s
    Chunksize: %s
""" % (inputfile, outputfolder, outfile_pattern, mode, UMIfile, logKept, logRejected, pairedEnd, cpuCount, chunk_size))

# Split the original file into chunks
bam = pysam.AlignmentFile(inputfile)
chunk = 0
totalReads = 0
reads_in_this_chunk = 0
old_name = None
outfile = pysam.AlignmentFile(tmp_pattern % chunk, "w", template = bam)

for read in bam.fetch(until_eof=True):
    totalReads += 1
    if old_name != read.query_name and reads_in_this_chunk > chunk_size:

        reads_in_this_chunk = 0
        chunk += 1
        outfile.close()
        outfile = pysam.AlignmentFile(tmp_pattern % chunk, "w", template = bam)
    outfile.write(read)
    old_name = read.query_name
    reads_in_this_chunk += 1
outfile.close()
bam.close()

number_of_chunks = chunk

# Make a list of all temporary files
tmp_files = []
for i in range(0, number_of_chunks+1):
    tmp_files.append(tmp_pattern % i)

# Extract UMIs from all chunks
def work(inBAM, outBAM):
    global chunk
    global listKept
    global listRejected

    reads_in_this_chunk = 0
    bam = pysam.AlignmentFile(inBAM)
    outfile = pysam.AlignmentFile(outBAM, "w", template = bam)
    for read in bam.fetch(until_eof=True):
        if reads_in_this_chunk % 100000 == 0:
            print("%s lines processed on %s ..." % (reads_in_this_chunk, multiprocessing.current_process()))
        UMI_from_read = str(read).split("\t")[0].split("_")[1]
        if mode == "fixed":
            if str(UMI_from_read) in tuple(UMI_list_6nt[0]):
                outfile.write(read)
                listKept.append(UMI_from_read)
            else:
                listRejected.append(UMI_from_read)
        else:
            if "N" not in str(UMI_from_read):
                outfile.write(read)
                listKept.append(UMI_from_read)
            else:
                listRejected.append(UMI_from_read)
        old_name = read.query_name
        reads_in_this_chunk += 1
    bam.close()
    outfile.close()

# Process all chunks
with Manager() as manager:
    # Define empty lists
    chunk = 0
    old_name = None
    listKept = manager.list()
    listRejected = manager.list()
    pool=Pool(processes=cpuCount)
    bam_files_list = []
    for file in tmp_files:
        tmp_pattern2 = outputfolder + "/" + file.split('/')[-1].split('.')[0] + "_noNoise_tmp%d.bam" % chunk
        bam_files_list.append(tmp_pattern2)
        pool.apply_async(work, (file, tmp_pattern2))
        chunk += 1
    pool.close()
    pool.join()

    # Merge with samtools
    #bam_files_list = glob.glob(os.path.join(outputfolder, "*noNoise_tmp*.bam"))
    arg = "samtools merge -f %s %s" % (outfile_pattern, ' '.join(bam_files_list))
    arg = arg.split()
    subprocess.call(arg,shell=False)

    # Remove tmp files
    for f in bam_files_list:
        os.remove(f)
    for f in tmp_files:
        os.remove(f)

    end = time.time()
    time = end-start
    print(" ")
    print("Finished running. Runtime %.2f seconds." % time)

    dfKept = pd.DataFrame.from_dict(dict(Counter(listKept)), columns = ["counts"], orient = "index").sort_values(by = "counts", ascending = False)
    dfKept.index.name = "UMI"
    dfKept.to_csv(logKept)
    dfRejected = pd.DataFrame.from_dict(dict(Counter(listRejected)), columns = ["counts"], orient = "index").sort_values(by = "counts", ascending = False)
    dfRejected.index.name = "UMI"
    dfRejected.to_csv(logRejected)

    totalReadsConsidered = dfKept["counts"].sum() + dfRejected["counts"].sum()
    readsKept = dfKept["counts"].sum()
    readsRejected = dfRejected["counts"].sum()
    print("Total number of reads in bam file: %i" % (totalReads))
    print("Total number of reads considered: %i (%.2f%%)" % (totalReadsConsidered, 100*(totalReadsConsidered/totalReads)))
    print("Total number of reads kept: %i (%.2f%%)" % ((readsKept), 100*(readsKept/(totalReads))))
    print("Total number of reads rejected: %i (%.2f%%)" % ((readsRejected), 100*(readsRejected/(totalReads))))
