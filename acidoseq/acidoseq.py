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

from termcolor import colored
from colorama import init # use Colorama to make Termcolor work on Windows too

import click

###############################################################################

def load_taxondump(idpath):
    """Importing the Acidobacteria taxon IDs"""
    taxons = {}
    with open(idpath) as csvfile:
        reader = csv.reader(csvfile, delimiter=',') 
        for row in reader:
            taxons[row[1]] = row[0]
        return taxons

def taxon_file(taxdumptype):
    while "ALL" not in taxdumptype and "U" not in taxdumptype:
        print("Error...")
        taxdumptype = input("Input here ('ALL' or 'U')?: ")
    if taxdumptype == "ALL":
        idpath = "input/acido_taxid_all.csv"
        taxons = load_taxondump(idpath) 
    if taxdumptype == "U":
        idpath = "input/acido_taxid_unclassified.csv"
        taxons = load_taxondump(idpath) 
    return taxons

###############################################################################
# Kaiju input and Acidobacteria read outputs

def insert_csv(kaijufile):
    """Inserting the CSV file of your results and returning a dictionary of them."""
    dict_seqid_taxon = {}
    with open(kaijufile) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            if row[1] not in dict_seqid_taxon:
                dict_seqid_taxon[row[1]] = { "reads": [] }     
            dict_seqid_taxon[row[1]]["reads"].append(row[0])
    return dict_seqid_taxon # above swaps over the columns and avoids the duplicates

def percentage(part, whole):
    """Calculating the coverage of Acidobacteria reads from the set of sequences."""
    return 100 * float(part)/float(whole)

###############################################################################
# ACGT comparisons

def list_of_sequences(fasta): # list of sequences, e.g. ['ACACCT', 'TAGC']
    reads = []
    for read in fasta.references:
        reads.append(fasta.fetch(read))
    return reads


def calculate_at(read):
    """Returns AT count."""
    return (read.lower().count("a") + read.lower().count("t")) / len(read) * 100.0
def calculate_gc(read):
    """Returns GC count."""
    return (read.lower().count("g") + read.lower().count("c")) / len(read) * 100.0

def plot_hist(myDict, style, taxdumptype): 
    plt.style.use(style)
    #max_read = int(max(myDict.values()))
    #min_read = int(min(myDict.values())) # random.randint(min_read, max_read)
        
    plt.hist(myDict.values(), bins=1000) 
    plt.xlabel('Ratio')   
    plt.ylabel('Count')
    
    meanDict = (sum(myDict.values())/float(len(myDict.values())))
    plt.axvline(x=meanDict, color='k')
    plt.text(x=(meanDict+1), y=(random.randint(20, 200)), s=str("%.2f" % meanDict))

    if taxdumptype == "U":
        ttype = "unclassified"
        plt.title('Histogram of ACGT for a collection of\n%s Acidobacteria sequences' % (ttype))
    elif taxdumptype == "ALL":
        ttype = "all"
        plt.title('Histogram of ACGT for a collection of\nAcidobacteria sequences')

    plt.grid(True)
    plt.savefig('output/acgt-comparison_%s_style-%s_%s.png' % (ttype, str(style), time_stamp))

###############################################################################
# GC ratio

def plot_hist_gc(myDict, style, ph, plottype, taxdumptype): 
    """Returns a plot of the GC ratio for a series of Acidobacteria sequences. 
    Includes the averages of the subdivisions based on the pH number."""
    plt.style.use(style) 
    ph = float(ph)
    
    plt.hist(myDict.values(), bins=1000, color="grey") 
    plt.xlabel('GC Ratio')   
    plt.ylabel('Count')  
    
    if plottype == "line": # mean GC
        lowph = {
            "sub1":58,
            "sub2":57.5,
            "sub3":59,
            "sub13":58.5
        }
        highph = {
            "sub4":60.12,
            "sub6":67.5,
            "sub22":67
        }
        medph = {
            "sub5":65.43, # medium
            "sub8":66.84, # medium
            "sub23":63 # medium
        } 
        colours = ['b', 'r', 'y', 'g']
        if ph < 5:
            for si, sub in enumerate(lowph):
                plt.axvline(lowph[sub], color=colours[si])
                plt.text(x=(lowph[sub]+0.5), y=(random.randint(20, 200)), s=str(sub), color=colours[si])
        elif ph > 5:
            for si, sub in enumerate(highph):
                plt.axvline(highph[sub], color=colours[si])
                plt.text(x=(highph[sub]+0.5), y=(random.randint(20, 200)), s=str(sub), color=colours[si])
        elif ph == 5:
            for si, sub in enumerate(medph):
                plt.axvline(medph[sub], color=colours[si])
                plt.text(x=(medph[sub]+0.5), y=(random.randint(20, 200)), s=str(sub), color=colours[si])    
                
    elif plottype == "span": # span of GC 
        if ph < 5:
            plt.axvspan(35.18, 67.1, alpha=0.25, color='green') # 1 
            plt.text(x=(35.18), y=(random.randint(20, 200)), s="sub1", color='green')
            plt.axvspan(57, 58, alpha=0.25, color='yellow') # 2 
            plt.text(x=(57.5), y=(random.randint(20, 200)), s="sub2", color='yellow')
            plt.axvspan(51, 73, alpha=0.25, color='red') # 3
            plt.text(x=(51), y=(random.randint(20, 200)), s="sub3", color='red')
            plt.axvspan(58, 59, alpha=0.25, color='blue') # 13
            plt.text(x=(58), y=(random.randint(20, 200)), s="sub13", color='blue')
        elif ph > 5:
            plt.axvspan(50, 61, alpha=0.25, color='blue') # 4 
            plt.text(x=(50), y=(random.randint(20, 200)), s="sub4", color='blue')
            plt.axvspan(67, 68, alpha=0.25, color='darkgreen') # 6 
            plt.text(x=(67), y=(random.randint(20, 200)), s="sub6", color='green')
            plt.axvspan(66, 67.5, alpha=0.25, color='pink') # 22 
            plt.text(x=(66), y=(random.randint(20, 200)), s="sub22", color='red')
        elif ph == 5:
            plt.axvspan(62.3, 68.3, alpha=0.25, color='red') # 5
            plt.text(x=(62.5), y=(random.randint(20, 200)), s="sub5", color='red')
            plt.axvspan(55.14, 71.83, alpha=0.25, color='pink') # 8
            plt.text(x=(55.5), y=(random.randint(20, 200)), s="sub8", color='pink')
            plt.axvspan(62, 64, alpha=0.25, color='orange') # 23        
            plt.text(x=(62), y=(random.randint(20, 200)), s="sub23", color='orange')
    if taxdumptype == "U":
        ttype = "unclassified"
        plt.title('Histogram of GC ratio of pH%.2f for a\ncollection of %s Acidobacteria sequences' % (ph, ttype))
    elif taxdumptype == "ALL":
        ttype = "all"
        plt.title('Histogram of GC ratio of pH%.2f for a\ncollection of Acidobacteria sequences' % (ph))
    plt.savefig('output/gc-ratio_%s_ph%.2f_plot-%s_style-%s_%s.png' % (ttype, ph, plottype, style, time_stamp))

###############################################################################

@click.command()
@click.option('--taxdumptype', default='ALL', help='Study "ALL" or only unclassified "U"?')
@click.option('--kaijufile', help='Place edited Kaiju (csv) in directory for ease.')
@click.option('--fastapath', help='Place FASTA in directory for ease.')
@click.option('--style', help="['seaborn-bright', 'seaborn-poster', 'seaborn-white', 'bmh', 'seaborn-darkgrid', 'seaborn-pastel', 'grayscale', '_classic_test', 'ggplot', 'seaborn-whitegrid', 'seaborn-dark', 'seaborn-muted', 'seaborn-colorblind', 'seaborn-ticks', 'Solarize_Light2', 'seaborn-notebook', 'dark_background', 'fast', 'seaborn', 'fivethirtyeight', 'seaborn-paper', 'seaborn-dark-palette', 'seaborn-talk', 'classic', 'seaborn-deep']")
@click.option('--plottype', default='line', help='"span" range of GC means OR "line" average mean GC')
@click.option('--ph', default='5', help='pH of soil, use map script for assistance.')
def main(taxdumptype, kaijufile, fastapath, style, plottype, ph):
    colors = ["yellow", "magenta", "blue", "green"]
    ph = float(ph)
    
    taxons = taxon_file(taxdumptype) # ALL/U
    taxon_read_map = insert_csv(kaijufile) # KAIJU
    
    has_taxon = 0
    total_reads = 0
    numrec = 0
    acido_reads = []
    
    for taxon_id in taxon_read_map:
        numrec += 1
        total_reads += len(taxon_read_map[taxon_id]["reads"])
        print("Record\t%s" % (str(numrec)))
        try:
            taxon_read_map[taxon_id]["scientific_name"] = taxons[taxon_id]
            has_taxon += len(taxon_read_map[taxon_id]["reads"])
            acido_reads.extend(taxon_read_map[taxon_id]["reads"])
        except KeyError:
            continue
        
    print(colored('\nAcidobacteria coverage of file:', colors[0])) 
    acido_coverage = percentage(has_taxon, total_reads)
    print("%.2f%%" % (acido_coverage))

############################################
    #FASTA
    fasta = pysam.FastaFile(fastapath) # FASTA
    output_acido_file = ("output/acido-reads_%s_%s.fa" % (taxdumptype, time_stamp))
    with open(output_acido_file, "w") as output:
        for r in acido_reads:                                                                                                
            seq = fasta.fetch(reference=r) 
            output.write(">%s\n%s\n" % (r, seq))
    print(colored("\nSuccessful! You can find the file here:", colors[1])) 
    print("%s" % (output_acido_file))
    
############################################
    # ACGT comparisons 
    
    out_path = output_acido_file
    fasta = pysam.FastaFile(out_path) 

    reads = list_of_sequences(fasta)   
    #max_read = len(max(reads, key=len))
    #min_read = len(min(reads, key=len))
    lens = [len(x) for x in reads]
    max_read = max(lens)
    min_read = min(lens)
    
    print(colored("\nStatistics:", colors[2]))
    print("Read Lengths\tMin: %d\tMax: %d" % (min_read, max_read))
    
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

############################################

    styles = plt.style.available
    try:
        plt.style.use(style)
    except OSError:
        print("\nStyle unrecognised, automatically choosing one for you...")
        style_int = random.randint(0,24)
        style = styles[style_int]
    
    x = 1
    plt.figure(x)
    myDict = at    
    plot_hist(at, style, taxdumptype)
    myDict = gc
    plot_hist(gc, style, taxdumptype)
    
############################################
    # GC ratio
    
    x = 2
    plt.figure(x)    
    plot_hist_gc(gc, style, ph, plottype, taxdumptype)
    
    
    print(colored("\nAll Done!\n", colors[3]))
    
###############################################################################
###############################################################################
###############################################################################

if __name__ == "__main__":
    time_stamp = strftime("%Y-%m-%d_%H-%M-%S", gmtime())
    main()    
