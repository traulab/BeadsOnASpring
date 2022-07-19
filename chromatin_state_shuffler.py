# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 14:04:05 2019

@author: s4138855
"""

from random import shuffle

D_compartments = {}
D_regions = {}

with open("wgEncodeBroadHmmGm12878HMM.bed","r") as F_chromatin:
    for line in F_chromatin:
        line = line.strip().split("\t")
        
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        comp = line[3]
        
        if comp == "NA":
            continue
        
        size = int(end - start)

        count = int(size/200)
        
        for x in range(0,count):
            if chrom in D_regions:
                D_regions[chrom].append((chrom,start+200*x))
            else:
                D_regions[chrom] = [(chrom,start+200*x)]

        if chrom in D_compartments:
            D_compartments[chrom].append((size, comp))
        else:
            D_compartments[chrom] = [(size, comp)]



with open("wgEncodeBroadHmmGm12878HMM_shuffled.bed","w") as F_chromatin:
    for chrom in D_regions:
        shuffle(D_compartments[chrom])  
        prev_end = D_regions[chrom][0][1]

        for i, region in enumerate(D_regions[chrom]):
            if region[1] >= prev_end:
                try:
                    F_chromatin.write(region[0] + "\t" + str(region[1]) + "\t" 
                                      + str(region[1] + D_compartments[chrom][0][0]) + "\t" 
                                      + D_compartments[chrom][0][1] + "\n")
                except:
                    break
                
                prev_end = region[1] + D_compartments[chrom][0][0]
                D_compartments[chrom].pop(0)
                    
                    
                    
