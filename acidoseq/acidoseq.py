# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 11:56:11 2018
@author: samantha
"""

import csv                                                                                                        
import pysam  
from time import gmtime, strftime
import collections
import matplotlib.pyplot as plt
import random

# Kaiju input and Acidobacteria read outputs

def insert_csv(file):
    """Inserting the CSV file of your results and returning a dictionary of them."""
    dict_seqid_taxon = {}
    with open(file) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            if row[1] not in dict_seqid_taxon:
                dict_seqid_taxon[row[1]] = { "reads": [] }     
            dict_seqid_taxon[row[1]]["reads"].append(row[0])
    return dict_seqid_taxon # above swaps over the columns and avoids the duplicates

def load_taxondump(path):
    """Importing the Acidobacteria taxon IDs"""
    taxons = {}
    with open(path) as csvfile:
        reader = csv.reader(csvfile, delimiter=',') 
        for row in reader:
            taxons[row[1]] = row[0]
        return taxons

def percentage(part, whole):
    """Calculating the coverage of Acidobacteria reads from the set of sequences."""
    return 100 * float(part)/float(whole)

######################################################

# ACGT comparisons

def list_of_sequences(fasta): # list of sequences, e.g. ['ACACCT', 'TAGC']
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

def calculate_at(seq):
    """Returns AT count."""
    return (seq.lower().count("a") + seq.lower().count("t")) / len(seq) * 100.0
def calculate_gc(seq):
    """Returns GC count."""
    return (seq.lower().count("g") + seq.lower().count("c")) / len(seq) * 100.0

def plot_hist(myDict, style): 
    plt.style.use(style)
    plt.hist(myDict.values(), bins=1000) 
    plt.xlabel('Ratio')   
    plt.ylabel('Count')
    
    meanDict = (sum(myDict.values())/float(len(myDict.values())))
    plt.axvline(x=meanDict, color='k')
    plt.text(x=(meanDict+1), y=(random.randint(50,500)), s=str("%.2f" % meanDict))
    
    plt.title('Histogram of ACGT for a collection of\nAcidobacteria sequences')
    plt.grid(True)
    #plt.show()
    plt.savefig('output/acgt-comparison_style-%s_%s.png' % (style, time_stamp))

######################################################

# GC Ratio

def plot_hist_gc(myDict, style, ph, plot_type, tax_type): 
    """Returns a plot of the GC ratio for a series of Acidobacteria sequences. 
    Includes the averages of the subdivisions based on the pH number."""
    ph = float(ph)
    plt.style.use(style)
    plt.hist(myDict.values(), bins=1000, color="grey") 
    plt.xlabel('GC Ratio')   
    plt.ylabel('Count')  
    
    if plot_type == "line":
        lowph = {
            "sub1":58,
            "sub2":57.5,
            "sub3":62,
            "sub13":58.5
        }
        highph = {
            "sub4":60,
            "sub6":68.5,
            "sub22":67
        }
        medph = {
            "sub5":65.5, # medium
            "sub8":66.8, # medium
            "sub23":63 # medium
        } 
        colours = ['b', 'r', 'y', 'g']
        if ph < 5:
            for si, sub in enumerate(lowph):
                plt.axvline(lowph[sub], color=colours[si])
                plt.text(x=(lowph[sub]+0.5), y=(random.randint(50,500)), s=str(sub))
        elif ph > 5:
            for si, sub in enumerate(highph):
                plt.axvline(highph[sub], color=colours[si])
                plt.text(x=(highph[sub]+0.5), y=(random.randint(50,500)), s=str(sub))
        elif ph == 5:
            for si, sub in enumerate(medph):
                plt.axvline(medph[sub], color=colours[si])
                plt.text(x=(medph[sub]+0.5), y=(random.randint(50,500)), s=str(sub))    
    elif plot_type == "span":
        if ph < 5:
            plt.axvspan(35.18, 57, alpha=0.5, color='green') # 1 
            plt.text(x=(35.18), y=(random.randint(50,500)), s="sub1")
            plt.axvspan(57.5, 57.7, alpha=0.5, color='yellow') # 2 
            plt.text(x=(57.5), y=(random.randint(50,500)), s="sub2")
            plt.axvspan(61, 63, alpha=0.5, color='blue') # 3
            plt.text(x=(61), y=(random.randint(50,500)), s="sub3")
            plt.axvspan(58, 59, alpha=0.5, color='darkblue') # 13
            plt.text(x=(58), y=(random.randint(50,500)), s="sub13")
        elif ph > 5:
            plt.axvspan(50, 60, alpha=0.5, color='yellow') # 4 
            plt.text(x=(50), y=(random.randint(50,500)), s="sub4")
            plt.axvspan(67, 68, alpha=0.5, color='green') # 6 
            plt.text(x=(67), y=(random.randint(50,500)), s="sub6")
            plt.axvspan(66, 67.5, alpha=0.5, color='darkgreen') # 22 
            plt.text(x=(66), y=(random.randint(50,500)), s="sub22")
        elif ph == 5:
            plt.axvspan(62.5, 68, alpha=0.5, color='orange') # 5
            plt.text(x=(62.5), y=(random.randint(50,500)), s="sub5")
            plt.axvspan(55, 66, alpha=0.5, color='pink') # 8
            plt.text(x=(55), y=(random.randint(50,500)), s="sub8")
            plt.axvspan(62, 64, alpha=0.5, color='red') # 23        
            plt.text(x=(62), y=(random.randint(50,500)), s="sub12")
        
    if tax_type == "U":
        ttype = "unclassified"
        plt.title('Histogram of GC ratio of pH%.2f for a\ncollection of %s Acidobacteria sequences' % (ph, ttype))
    elif tax_type == "all":
        ttype = "all"
        plt.title('Histogram of GC ratio of pH%.2f for a\ncollection of Acidobacteria sequences' % (ph))
    #plt.grid(True)
    #plt.show()
    plt.savefig('output/gc-ratio_%s_ph%.2f_plot-%s_style-%s_%s.png' % (ttype, ph, plot_type, style, time_stamp))


####################################################################################################
####################################################################################################

if __name__ == "__main__":
    time_stamp = strftime("%Y-%m-%d_%H-%M-%S", gmtime()) 
    
    # Kaiju input and Acidobacteria read outputs    
    
    taxdump_type = input("All species or only unclassified ('all' or 'U')? ")
    if taxdump_type == "all":
        taxons = load_taxondump("input/acido_taxid_all.csv") 
    if taxdump_type == "U":
        taxons = load_taxondump("input/acido_taxid_unclassified.csv")  
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
        
    acido_coverage = percentage(has_taxon, total_reads)
    print("Acidobacteria coverage of file: %.2f" % (acido_coverage))
      
    all_fasta_path = input("Enter the FASTA file of all reads: ") # e.g. all_reads.fa
    fasta = pysam.FastaFile(all_fasta_path)
    output_acido_file = ("output/acido-reads_%s_%s.fa" % (taxdump_type, time_stamp))
    with open(output_acido_file, "w") as output:
        for r in acido_reads:                                                                                                   
            seq = fasta.fetch(reference=r) 
            output.write(">%s\n%s\n" % (r, seq))

    print(output_acido_file)
    
######################################################
    
    # ACGT comparisons 
    
    path = output_acido_file
    fasta = pysam.FastaFile(path) 

    reads = list_of_sequences(fasta)   
    max_read = len(max(reads, key=len))
    min_read = len(min(reads, key=len))
    print("Read Lengths\tMin: %i\tMax: %i" % (min_read, max_read))
    
    at = {}
    gc = {}
    for read in fasta.references:
        at[read] = calculate_at(fasta.fetch(read))  
        gc[read] = calculate_gc(fasta.fetch(read))
    max_at = max(at.values())
    min_at = min(at.values())  
    mean_at = (sum(at.values())/float(len(at.values())))
    print("AT\tMin: %f\tMax: %f\tMean: %f" % (min_at, max_at, mean_at))
    max_gc = max(gc.values())
    min_gc = min(gc.values())  
    mean_gc = (sum(gc.values())/float(len(gc.values())))
    print("GC\tMin: %f\tMax: %f\tMean: %f" % (min_gc, max_gc, mean_gc))

    print(plt.style.available)
    style = input("Insert the style you want: ")
    
    x = 1
    plt.figure(x)    
    plot_hist(at, style)
    plot_hist(gc, style)
    
######################################################
    
    # GC Ratio

    plot_type = input("Enter plot type ('span' or 'line'): ")
    ph = input("Insert pH of soil: ")
    x = 2
    plt.figure(x)    
    plot_hist_gc(gc, style, ph, plot_type, taxdump_type)