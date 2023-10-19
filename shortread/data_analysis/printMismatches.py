from Bio import pairwise2
import pandas as pd


UMI_file = pd.read_csv("/Users/rmvpaeme/test/UMI_trim/UMI_list_RNA_exome.txt", header = None)

UMI_list = list(UMI_file[0])

UMI_to_6nt = []
for UMI in seqs:
    if len(UMI) == 7:
        UMI = UMI[:6]
        UMI_to_6nt.append(UMI)
    elif len(UMI) == 6:
        UMI_to_6nt.append(UMI)

def levenshteinDistance(s1, s2):
    if len(s1) > len(s2):
        s1, s2 = s2, s1

    distances = range(len(s1) + 1)
    for i2, c2 in enumerate(s2):
        distances_ = [i2+1]
        for i1, c1 in enumerate(s1):
            if c1 == c2:
                distances_.append(distances[i1])
            else:
                distances_.append(1 + min((distances[i1], distances[i1 + 1], distances_[-1])))
        distances = distances_
    return distances[-1]


def hamming_distance(s1, s2):
    """Return the Hamming distance between equal-length sequences"""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))

for sequence in UMI_list:
    distances = []
    for sequence2 in UMI_list:
        if sequence2 != sequence:
            distances.append(levenshteinDistance(sequence2, sequence))
    print("Minimum levenshtein distance for %s: %i" % (sequence, min(distances)))
    print("Maximum levenshtein distance for %s: %i" % (sequence, max(distances)))

total_distances = []
for sequence in UMI_to_6nt:
    distances_per = []
    for sequence2 in UMI_to_6nt:
        if sequence2 != sequence:
            distances_per.append(hamming_distance(sequence2, sequence))
            total_distances.append(hamming_distance(sequence2, sequence))
    print("Minimum hamming distance for %s: %i" % (sequence, min(distances_per)))
    print("Maximum hamming distance for %s: %i" % (sequence, max(distances_per)))
print("Minimum total hamming distance for all 120 UMIs: %i" % min(total_distances))
print("Maximum total hamming distance for all 120 UMIs: %i" % max(total_distances))
