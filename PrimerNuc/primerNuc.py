# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 12:14:34 2018
@author: Jennifer Lu & Andrew Johnston
"""

from operator import itemgetter, is_not
from collections import OrderedDict
import datetime
import time
import sys
from math import log10, sqrt
import re
from functools import partial 
from time import gmtime, strftime
from importlib import import_module
import os
import argparse

#############################################################
def readBed2Dict(inputBedFile):
    
    BedDict = OrderedDict()
    with open(inputBedFile, 'r') as bedFile:
        for line in bedFile:
            #only reads lines which starts with chr
            if line.startswith("chr"):
                line = line.split()
    
                chrom = line[0]
                start = int(line[1])
                end = int(line[2])
                
                #check name of region
                try:
                    name = line[3]
                except:
                    name = ''.join([chrom, ':', line[1], '-', line[2]])
                
                rline = [chrom, start, end, name]
                #add to dict with the chrom as the key
                if chrom in BedDict:
                    BedDict[chrom].append(rline)
                else:
                    BedDict[chrom] = [rline]
            
    return BedDict

#############################################################
def Result2Report(NucResultTable, resultFile):

    with open(resultFile, 'w') as rfile:
    
        for line in NucResultTable:
            line = '\t'.join(line)+'\n'
            rfile.write(line)

    return rfile

#############################################################

#main function for nucleosome peaks
    
def Nucleosome(inputBedFile, resultFile):
    
   
    title = ['chrom', 'Peak', 'name', 'region', 'WPS', "to 5'end", "to 3' end"]
    
    resultsTable = [title]

    BedDict = readBed2Dict(inputBedFile) #dictionary of input nucleosome regions

    for chrom in BedDict:
        regions = BedDict[chrom]
        
        if chrom == 'chr1':
            from Chrom1peaks import getPeaks
        
        elif chrom == 'chr2':
            from Chrom2peaks import getPeaks
        
        elif chrom == 'chr3':
            from Chrom3peaks import getPeaks
            
        elif chrom == 'chr4':
            from Chrom4peaks import getPeaks
            
        elif chrom == 'chr5':
            from Chrom5peaks import getPeaks
        
        elif chrom == 'chr6':
            from Chrom6peaks import getPeaks
            
        elif chrom == 'chr7':
            from Chrom7peaks import getPeaks
            
        elif chrom == 'chr8':
            from Chrom8peaks import getPeaks
            
        elif chrom == 'chr9':
            from Chrom9peaks import getPeaks
            
        elif chrom == 'chr10':
            from Chrom10peaks import getPeaks
            
        elif chrom == 'chr11':
            from Chrom11peaks import getPeaks
        
        elif chrom == 'chr12':
            from Chrom12peaks import getPeaks
        
        elif chrom == 'chr13':
            from Chrom13peaks import getPeaks
            
        elif chrom == 'chr14':
            from Chrom14peaks import getPeaks
            
        elif chrom == 'chr15':
            from Chrom15peaks import getPeaks
        
        elif chrom == 'chr16':
            from Chrom16peaks import getPeaks
            
        elif chrom == 'chr17':
            from Chrom17peaks import getPeaks
            
        elif chrom == 'chr18':
            from Chrom18peaks import getPeaks
            
        elif chrom == 'chr19':
            from Chrom19peaks import getPeaks
            
        elif chrom == 'chr20':
            from Chrom20peaks import getPeaks

        elif chrom == 'chr21':
            from Chrom21peaks import getPeaks
            
        elif chrom == 'chr22':
            from Chrom22peaks import getPeaks
            
        elif chrom == 'chrX':
            from ChromXpeaks import getPeaks
            
        elif chrom == 'chrY':
            from ChromYpeaks import getPeaks
        else:
            continue


        ###
        #go through each line in the region
                
        if regions:
            for line in regions:
                chr = line[0]
                start = (line[1] + line[2])//2
                end = start + 1
                name = line[3]

                PeakDict =  getPeaks(start, end) #dict of peaks within the region

                value = []
                for coord in sorted(PeakDict):
                    rline = PeakDict[coord]
#                    print(rline)
                    region = rline[0]
                    score = rline[1]
                    to5p = rline[2]
                    to3p = rline[3]

                    if args.NPP_direction == 'up':
                        if name == '+' and region == "3' flanking":
                            value = [chr, coord, coord + 1,name, region, score, to5p, to3p]
                            value = [str(n) for n in value]
                        elif name == '-' and region == "5' flanking":
                            value = [chr, coord, coord + 1,name, region, score, to5p, to3p]
                            value = [str(n) for n in value]
                    elif args.NPP_direction == 'down':
                        if name == '+' and region == "5' flanking":
                            value = [chr, coord, coord + 1,name, region, score, to5p, to3p]
                            value = [str(n) for n in value]
                        elif name == '-' and region == "3' flanking":
                            value = [chr, coord, coord + 1,name, region, score, to5p, to3p]
                            value = [str(n) for n in value]
                    elif args.NPP_direction == 'nearest':
                        if value == []:
                            value = [chr, coord, coord + 1,name, region, score, to5p, to3p]
                            value = [str(n) for n in value]
                            prev_to5p = to5p
                        elif to5p < prev_to5p:
                            value = [chr, coord, coord + 1,name, region, score, to5p, to3p]
                            value = [str(n) for n in value]
                            prev_to5p = to5p
                    else:
                        print("direction value not assigned: down, up or nearest")
                        exit(1)

                resultsTable.append(value)
 
        #############
        #remove module get peaks after use
        del(getPeaks) 


    #write results to report
    result = Result2Report(resultsTable, resultFile)

    
    return resultsTable

parser = argparse.ArgumentParser()
parser.add_argument('output_bed', help='input bed file')
parser.add_argument('input_bed', help='output bed file')
parser.add_argument('NPP_direction', help='direction of nearest NPP (up,down,or either')
args = parser.parse_args()

Nucleosome(args.output_bed, args.input_bed)

