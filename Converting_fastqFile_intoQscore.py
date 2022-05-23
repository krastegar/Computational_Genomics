#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 17:43:27 2022

@author: bioinfo
"""

from Bio import SeqIO
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np 


qscore = []
for record in SeqIO.parse("assignment1.fastq", "fastq"): #using biopython to parse through entire fastq file
    scores = record.letter_annotations["phred_quality"] # variable that shows all of the phred quality scores
    qscore.append(scores)#appending scores to an array
    #print(scores)
df = pd.DataFrame(qscore)#transofrming array into dataframe
l = len(df.T)+1

#-------Customizing my plots--------------------------
f,ax=plt.subplots(figsize=(10,5)) #f, has the same role as plt.figure(figsize = ())
df.mean().plot(ax=ax,c='black') #here we are looking at a smooth curve going through the mean at each position
boxprops = dict(linestyle='-', linewidth=1, color='black') 
ax.set_xticks(np.arange(0, l, 5))
ax.set_xticklabels(np.arange(0, l, 5))
ax.set_xlabel('position(bp)')
ax.set_xlim((0,l))
ax.set_ylim((0,40))
ax.set_title('per base sequence quality')    

