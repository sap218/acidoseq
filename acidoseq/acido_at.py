# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 12:50:13 2018

@author: samantha
"""

import pysam
import collections
from time import gmtime, strftime

def calculate_at(seq): # ratio at AT in each sequence
    return (seq.lower().count("a") + seq.lower().count("t")) / len(seq) * 100.0

def list_of_sequences(fasta): # list of sequences, e.g. ['AGCT', 'TAGC']
    reads = []
    for read in fasta.references:
        reads.append(fasta.fetch(read))
    return reads

def counting_acgt(list): # getting ACGT count for each sequence
    results = []
    x = 0
    for letter in list:
        result = collections.Counter(list[x])
        x = x + 1
        results.append(result)
    return results

##########

if __name__ == "__main__":
    path = input("enter fasta file: ") # acido_reads_2018-07-28_22-28-17.fa
    fasta = pysam.FastaFile("/home/samantha/Dissertation/python/%s" % path) 

    reads = list_of_sequences(fasta)
    #results = counting_acgt(reads)

    at = {}
    for read in fasta.references:
        at[read] = calculate_at(fasta.fetch(read))
    for key, value in at.items():
        if value > 55:
            print("%s\t%.2f" % (key, value))
    '''   
    time_stamp = strftime("%Y-%m-%d_%H-%M-%S", gmtime())    
    with open("acido_high-at-reads_%s.fa" % (time_stamp), "w") as output:
        for r in at:                    
            output.write("%s\n" % (r))
    '''            
            