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
    return (seq.lower().count("g") + seq.lower().count("c")) / len(seq) * 100.0

def plot_hist(myDict, style): 
    subdivisions = {"sub1":58, "sub2":57.6, "sub3":62,"sub4":60,"sub5":65.5,"sub6":67.2,"sub8":66.8,"sub13":58.5,"sub22":67.15,"sub23":63}
    
    plt.style.use(style)
    plt.hist(myDict.values(), bins=1000) 
    plt.xlabel('GC Ratio')   
    plt.ylabel('Count')
    
    for sub in subdivisions:
        plt.axvline(subdivisions[sub], color='k')
        plt.text(x=(subdivisions[sub]), y=(random.randint(100,1300)), s=str(sub))
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
        
    plt.title('Histogram of GC ratio for a collection\nof Acidobacteria sequences')
    plt.grid(True)
    #plt.show()
    plt.savefig('gc-ratio_%s.png' % time_stamp)
    
#######
    
if __name__ == "__main__":
    time_stamp = strftime("%Y-%m-%d_%H-%M-%S", gmtime()) 
    path = input("Enter FASTA file: ") # e.g. acido_reads_2018-07-28_22-28-17.fa
    fasta = pysam.FastaFile(path) 
    
    gc = {}
    for read in fasta.references:
        gc[read] = calculate_gc(fasta.fetch(read))   
        
    max_gc = max(gc.values())
    min_gc = min(gc.values())  
    mean_gc = (sum(gc.values())/float(len(gc.values())))
    print("Min\t%f\nMax\t%f\nMean\t%f" % (min_gc, max_gc, mean_gc))

    print(plt.style.available)
    style = input("insert style you want: ")
    x = 1
    plt.figure(x)    
    plot_hist(gc, style)
    
    #########
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