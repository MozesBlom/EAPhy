#Required modules
import os
from Bio import SeqIO

def getContigStats(locus, locus_path):
    bestContig_file = open (locus_path, 'r')
    bestContig = SeqIO.parse(bestContig_file, 'fasta')
    indiv_count = 0
    indiv_length = []
    total_length = 0
    gap = "-"
    for record in bestContig:
        indiv_count += 1
        indiv_length = len(record.seq) - (record.seq.count(gap, start = 0, end = len(record.seq)))
        total_length = total_length + indiv_length
    if indiv_count == 0:
        average_length = 0
    else:
        average_length = total_length / indiv_count
    return locus, indiv_count, average_length
