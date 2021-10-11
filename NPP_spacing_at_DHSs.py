"""
Created on Tue Dec 22 11:00:01 2020

@author: Andrew Douglas Johnston

This script plots the spacing of nucleosome protection peaks (NNPs) flanking the centers
of DHS sites of ENCODE cell lines and calculates the ratio of 250bp/185bp spacing

Inputs include:
    NPP positions in BED format converted from big beds (e.g. GSE71378_CH01.bed)
    DHS site positions in BED format converted from big beds (e.g. wgEncodeDukeDnaseGM12891.fdr01peaks.hg19.bb.bed)


"""

from matplotlib import pyplot
import numpy as np
import os

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


L_filenames = ['GSE71378_CH01.bed','GSE71378_BH01.bed','BBC01_NPPs.bed',
               'GSE71378_IH01.bed','BRE01_NPPs.bed','COL01_NPPs.bed',
               'BRA01_NPPs.bed','GSE71378_IH02.bed',
               'GSE71378_CA01.bed']

for NPP_filename in L_filenames:
    D_NPP_sites = {}
       
    with open(NPP_filename, 'r') as F_NPP_bed:
        for line in F_NPP_bed:
            NPP = line.strip().split("\t")
            chromosome = NPP[0]
            NPP_position = int(NPP[6])
    
            if chromosome not in D_NPP_sites:
                D_NPP_sites[chromosome] = [NPP_position]
            else:
                D_NPP_sites[chromosome].append(NPP_position)
    
    rootdir = os.path.dirname(os.path.realpath(__file__))
    
    F_ratios = open(NPP_filename + '_ratios.txt','w')
    F_ratios.close() 
    
    counter = 0
    for subdir, dirs, files in os.walk(rootdir):
    
        for file in files:
    
            #print os.path.join(subdir, file)
            filepath = subdir + os.sep + file
    
            if filepath.endswith(".bb.bed"):
                print (filepath)
                
                counter += 1
                
                # if counter > 1:
                    # break
    
                D_DHS_sites = {}
            
                with open(filepath, 'r') as F_DHS_bed:
                    for line in F_DHS_bed:
                        DHS = line.strip().split("\t")
                        chromosome = DHS[0]
                        start_nuc = int(DHS[1])
                        end_nuc = int(DHS[2])
                        
                        DHS_center = int((start_nuc + end_nuc)/2)
                        
                        if chromosome not in D_DHS_sites:
                            D_DHS_sites[chromosome] = [DHS_center]
                        else:
                            D_DHS_sites[chromosome].append(DHS_center)
            
            
                D_DHS_matched_NPP = {}
                            
                for chromosome in D_DHS_sites:
                    D_DHS_matched_NPP[chromosome] = []
                    i = 0
                    n = 1
                    while i < len(D_DHS_sites[chromosome]):
                        if n > len(D_NPP_sites[chromosome])-1:
                            break
                        
                        NPP = D_NPP_sites[chromosome][n-1]
                        NPP2 = D_NPP_sites[chromosome][n]
                        
                        if NPP <= D_DHS_sites[chromosome][i] <= NPP2:
                            D_DHS_matched_NPP[chromosome].append(NPP2 - NPP)
                            i += 1
                            n += 1
                        elif D_DHS_sites[chromosome][i] < NPP:
                            i += 1
                        else:
                            n += 1
                            
                
                with open(filepath +'_' + NPP_filename + '.txt', 'w') as F_output:   
                    L_nuc_distances = []
                    for chromosome in D_DHS_matched_NPP:
                        for x in D_DHS_matched_NPP[chromosome]:
                            F_output.write(str(x) + '\n')
                            L_nuc_distances.append(x)
                    
                    L_unique_distances  = list(set(L_nuc_distances))
                    L_unique_distances = [x for x in L_unique_distances if 0 <= x <= 600]
                    
                    L_distance_counts = []
                    for distance in L_unique_distances:
                        L_distance_counts.append(L_nuc_distances.count(distance))
                    
                    L_relative_freq = []
                    for count in L_distance_counts:
                        L_relative_freq.append(count/sum(L_distance_counts))
                    
                    xhat = smooth(L_relative_freq,30)
                    L_unique_distances, xhat = zip(*sorted(zip(L_unique_distances, xhat)))
                    
                    pyplot.plot(L_unique_distances, xhat, lw=0.5)
    
    
                    index_185 = L_unique_distances.index(185)
                    index_250 = L_unique_distances.index(250)
                    cellLine = filepath.split('Dnase')[1].split('.')[0]
                    
                    F_ratios = open(NPP_filename + '_ratios.txt','a')
                    F_ratios.write(cellLine + '\t' +
                                   str(xhat[index_185]) + '\t' +
                                   str(xhat[index_250]) + '\t' +
                                   str(xhat[index_250]/xhat[index_185]) + '\n')
                    F_ratios.close() 
                    
    fig = pyplot.gcf()
    fig.set_size_inches(6,4, forward = False)
    #fig.set_size_inches(2.5,2.5, forward = False)
                    
    pyplot.ylim((-0.0001,0.015))
    # pyplot.xticks(ticks = np.arange(-4.5,400,21))
    pyplot.xlim((0,500))
        
    pyplot.savefig(NPP_filename +'_' + NPP_filename + '.svg')
    pyplot.show()                    