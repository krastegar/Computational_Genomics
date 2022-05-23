#------------------------------------------------------------------------------------------------------------
# Purpose: To simulate genetic drift of a population with inbreeding. Create a visual representation 
#          of what happens to the homozyogosity of the population when inbreeding is present? 
#------------------------------------------------------------------------------------------------------------
# Algorithm: 
#   Step 1: Declare variables (i.e, pop size, # of generations, etc...)
#   Step 2: Create a for loop which we simulate the genetic drift for each genotype
#           with inbreeding 
#   Step 3: Append probabilities to list to be used for plotting 
#   Step 4: Update probablilites inside of loop 
#

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from decimal import *
import random
from math import comb

p = 0.2 # allele frequency
q = 1-p
f = 0.2 # inbreeding coefficient 
n = 1000 # population size    x = random.randint(0,2*n)
generation = 50
p_AA,p_AG,p_GG = [],[],[]
col = ['AA','GG','AG']

for i in range(0, generation):
    prob_nextGen_AA = pow(Decimal(p),2) + Decimal(p)*Decimal(q)*Decimal(f)
    prob_nextGen_AG = 2*Decimal(p)*Decimal(q)*(1-Decimal(f))
    prob_nextGen_GG = pow(Decimal(q),2) + Decimal(p)*Decimal(q)*Decimal(f)
    p_AA.append(prob_nextGen_AA)
    p_GG.append(prob_nextGen_GG)
    p_AG.append(prob_nextGen_AG)
    p = (prob_nextGen_AA+prob_nextGen_AG)/2 # genetic drift under WF process
    x = random.randint(0,2*n)
    q = 1-p

# Formatting my dataframe
genotype = np.array((p_AA,p_GG,p_AG))
genotype = np.transpose(genotype)
genotype_freq = pd.DataFrame(genotype, columns=col)
number_generation = range(0,generation)

#-------Customizing my plots--------------------------
plt.figure()
plt.bar(number_generation, genotype_freq.iloc[:, 2], color='g', align='edge')
plt.bar(number_generation, genotype_freq.iloc[:, 1] , bottom =genotype_freq.iloc[:,2], color='purple', align = 'edge')
plt.bar(number_generation, genotype_freq.iloc[:, 0], bottom=genotype_freq.iloc[:, 2], color='b', align = 'edge')
plt.legend(col,loc="upper left")
plt.xlabel("Generations")
plt.ylabel("Allele Frequency")
plt.show()
