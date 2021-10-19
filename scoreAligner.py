
"""
Created on Fri Apr  2 14:37:02 2021

@author: Andrew Johnston

This script uses bigwig files containing genome-wide Windowed Protection Scores (WPS)
and BED files containing transcription factor (TF) binding sites
to produce plots of average WPS against distance from TF binding sites

"""

import numpy as np
import pyBigWig as bw
import matplotlib.pyplot as plt

chromo_range_start = 1
chromo_range_end = 25

window_size = 5000

D_TF_sites = {}

with open("TF_binding_sites.bed") as TF_bed:
    for line in TF_bed:

        chunk = line.strip().split('\t')
        
        chromo = chunk[0]
        try:
            if chromo != 'chrX' and chromo != 'chrY':
                int(chromo.split('chr')[1])
        
            if chromo not in D_TF_sites:
                D_TF_sites[chromo] = []
            
            start = int(chunk[1])
            end = int(chunk[2])
            strand = chunk[3]
            
            centre = (start + end)//2
            
            D_TF_sites[chromo].append((centre, strand))
        except:
            continue

TF_count = 0    
TF_l_wps = np.zeros(window_size, dtype=float)

for chromo in range(chromo_range_start,chromo_range_end):
    l_wps = []
        
    if chromo == 23:
        chromo = 'X'
    elif chromo == 24:
        chromo = 'Y'
    
    bw_file_name = "SampleName_chr" + str(chromo) + "_l_wps.bw"
    
    bigwig = bw.open(bw_file_name)
    chromo_size = bigwig.header()['nBasesCovered']
 
    if chromo == 1:
        chromo_size = 249240753
    l_wps = bigwig.values(str(chromo), 0, chromo_size, numpy=True)   
    chromo = 'chr' + str(chromo) 
    
    if chromo not in D_TF_sites:
        bigwig.close()
        continue
    
    for TF_centre in D_TF_sites[chromo]:
        window_array = l_wps[TF_centre[0]-window_size//2:TF_centre[0]+window_size//2]
        
        if len(window_array) == window_size:
            TF_count += 1  #count to check against input TF BED file
            
            if TF_centre[1] == "+":
                TF_l_wps[0:window_size+1] += window_array
            else:
                TF_l_wps[0:window_size+1] += np.flipud(window_array)
                
    bigwig.close()

   
TF_l_wps[0:window_size+1] = TF_l_wps[0:window_size+1]/TF_count

fig, ax = plt.subplots(figsize=(22, 5))
x_coordinate = [ i -window_size//2 for i in range(len(TF_l_wps)) ]
plt.plot(x_coordinate,TF_l_wps)

plt.savefig('TF_adj_WPS_aggregates_all_chr.svg')
plt.show()

