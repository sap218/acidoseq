# -*- coding: utf-8 -*-
"""
Created on Sun Jul 22 17:41:24 2018

@author: samantha
"""

import pysam
import matplotlib.pylab as plt
from time import gmtime, strftime
import random

def calculate_gc(seq):
    """Returns GC count for a set of Sequences"""
    return (seq.lower().count("g") + seq.lower().count("c")) / len(seq) * 100.0

def plot_hist(myDict, style, ph): 
    """Returns a plot of the GC ratio for a series of Acidobacteria sequences. 
    Includes the averages of the subdivisions based on the pH number."""
    ph = float(ph)
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
    plt.style.use(style)
    plt.hist(myDict.values(), bins=1000, color="grey") 
    plt.xlabel('GC Ratio')   
    plt.ylabel('Count')   
    
    # this stuff is a working progress, please be patient :-)
    colours = ['b', 'r', 'y', 'g']
    if ph < 5:
        for si, sub in enumerate(lowph):
            plt.axvline(lowph[sub], color=colours[si])
            plt.text(x=(lowph[sub]+0.5), y=(random.randint(100,1300)), s=str(sub))
    elif ph > 5:
        for si, sub in enumerate(highph):
            plt.axvline(highph[sub], color=colours[si])
            plt.text(x=(highph[sub]+0.5), y=(random.randint(100,1300)), s=str(sub))
    elif ph == 5:
        for si, sub in enumerate(medph):
            plt.axvline(medph[sub], color=colours[si])
            plt.text(x=(medph[sub]+0.5), y=(random.randint(100,1300)), s=str(sub))    
    """
    plt.axvspan(53, 60, alpha=0.5, color='green') # 1
    plt.axvspan(57.5, 57.7, alpha=0.5, color='yellow') # 2
    plt.axvspan(51, 73, alpha=0.5, color='blue') # 3
    plt.axvspan(50, 52, alpha=0.5, color='yellow') # 4
    plt.axvspan(62, 67, alpha=0.5, color='orange') # 5
    plt.axvspan(66, 68, alpha=0.5, color='green') # 6
    plt.axvspan(55, 71, alpha=0.5, color='pink') # 8
    plt.axvspan(57, 59, alpha=0.5, color='darkblue') # 13
    plt.axvspan(66.8, 67.5, alpha=0.5, color='darkgreen') # 22
    plt.axvspan(62, 64, alpha=0.5, color='red') # 23        
    """
    plt.title('Histogram of GC ratio of pH%.2f for a\ncollection of Acidobacteria sequences' % (ph))
    #plt.grid(True)
    #plt.show()
    plt.savefig('gc-ratio_ph%.2f_style-%s__%s.png' % (ph, style, time_stamp))
    
#######
    
if __name__ == "__main__":
    time_stamp = strftime("%Y-%m-%d_%H-%M-%S", gmtime()) 
    path = input("Enter FASTA file: ") # e.g. acido_reads_2018-07-28_22-28-17.fa
    fasta = pysam.FastaFile(path) 
    
    gc = {}
    for read in fasta.references:
        gc[read] = calculate_gc(fasta.fetch(read))   

    ph = input("Insert pH of soil: ")
    print(plt.style.available)
    style = input("Insert style you want: ")
    x = 1
    plt.figure(x)    
    plot_hist(gc, style, ph)
    
    #########
    # Hopefully output txt of the sequences that reside in the subdivisons.
    '''
    dict_less_50 = {}
    dict_50_55 = {}
    dict_55_60 = {}
    dict_60_65 = {}
    dict_more_65 = {}
    
    for key,val in gc.items():
        if val < 50: 
            dict_less_50[key] = val
        elif val < 55 and val > 50:
            dict_50_55[key] = val
        elif val < 60 and val > 55:
            dict_55_60[key] = val
        elif val < 65 and val > 60:
            dict_60_65[key] = val
        elif val > 65:
            dict_more_65[key] = val
    '''