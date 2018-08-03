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
    acido_genomes = input("Enter Acidobacteria Genomes: ") # path to acido_taxid.csv e.g. /input/acido_taxid.csv
    taxons = load_taxondump(acido_genomes)     
    path = input("Enter your Kaiju Output (edited) file: ") # e.g. result_seqid_taxon.csv
    taxon_read_map = insert_csv(path)

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
    all_fasta_path = input("Enter the FASTA file of all reads: ") # e.g. all_reads.fa
    fasta = pysam.FastaFile(all_fasta_path)
    with open("acido_reads_%s.fa" % (time_stamp), "w") as output:
        for r in acido_reads:                                                                                                   
            seq = fasta.fetch(reference=r) 
            output.write(">%s\n%s\n" % (r, seq))
