# -*- coding: utf-8 -*-
"""
Created on Sat Jul 21 14:19:49 2018
@author: samantha
"""

import csv                                                                                                        
import pysam  
from time import gmtime, strftime

def insert_csv(file):
    dict_seqid_taxon = {}
    with open(file) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            if row[1] not in dict_seqid_taxon:
                dict_seqid_taxon[row[1]] = { "reads": [] }     
            dict_seqid_taxon[row[1]]["reads"].append(row[0])
    return dict_seqid_taxon # above swaps over the columns and avoids the duplicates

def load_taxondump(path):
    taxons = {}
    with open(path) as csvfile:
        reader = csv.reader(csvfile, delimiter=',') 
        for row in reader:
            taxons[row[1]] = row[0]
        return taxons

def percentage(part, whole):
    return 100 * float(part)/float(whole)

#####

if __name__ == "__main__":
    taxons = load_taxondump("/home/samantha/Dissertation/python/acido_taxid.csv")
    taxon_read_map = insert_csv("/home/samantha/Dissertation/python/result_seqid_taxon.csv")       
    
    has_taxon = 0
    total_reads = 0
    numrec = 0
    acido_reads = []
    
    for taxon_id in taxon_read_map:
        numrec += 1
        total_reads += len(taxon_read_map[taxon_id]["reads"])
        print("record", numrec)

        try:
            taxon_read_map[taxon_id]["scientific_name"] = taxons[taxon_id] # try to find taxon id for acidobacteria
            has_taxon += len(taxon_read_map[taxon_id]["reads"])
            acido_reads.extend(taxon_read_map[taxon_id]["reads"])
        except KeyError:
            continue
    print("acido coverage of file", percentage(has_taxon, total_reads))
    
    time_stamp = strftime("%Y-%m-%d_%H-%M-%S", gmtime())    
    fasta = pysam.FastaFile("/home/samantha/Dissertation/all.fa")
    with open("acido_reads_%s.fa" % (time_stamp), "w") as output:
        for r in acido_reads:                                                                                                   
            seq = fasta.fetch(reference=r) 
            output.write(">%s\n%s\n" % (r, seq))
