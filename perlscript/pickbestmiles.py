
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec as gs
import subprocess as sp
from scipy.stats import gaussian_kde

#from multiprocessing import Pool
#from multiprocessing import Process,Queues

# MILES=5, COELHO=4
LIB = 5

FINAL_OUTPUT = 'Avg_par.dat'

if (LIB==4):
        FINAL_OUTPUT = 'Avg_par_Coelho.dat'

if (LIB==5):
        FINAL_OUTPUT = 'Avg_par_Miles.dat' 
        
####################################################################
        
        
best = {}        

# from Miles Create miles dictionary 
libstar  = np.genfromtxt('names_miles91.lis', delimiter=" ", usecols=0, dtype=str) 
id_miles = np.genfromtxt('names_miles91.lis', delimiter=" ", usecols=1, dtype=int)  
NamesMiles = dict(zip(libstar, id_miles ))

    #open atmosphericParameters and build a best dictionary     
with open('atmosphericParameters.dat') as apar:
   lines = [raw.strip().split() for raw in apar]
   for n, cols in enumerate(lines):
      if len(cols) == 1:
          best[cols[0]] = {'LIB': '', 'id': '', 'rv': ''}
          best[cols[0]]['LIB'] = lines[n+1][0]
          
          #compera em input a id in best
for star in best.keys():
    if best[star]['LIB'] in NamesMiles :
        best[star]['id'] = NamesMiles[best[star]['LIB']]
    
    
    
    #create a templatesForCC.txt
with open('./templatesForCC.txt', 'w') as tpCC:
     tpCC.write(str(len(best.keys()))+'\n')
     n=1
     for star in best.keys():
         tpCC.write(star+' '+str(n)+' '+best[star]['LIB']+' '+str(best[star]['id'])+'\n')
         n = n + 1
         
    
    
# outputs templates for CC, 'best' is a dictionary contains best star, ( the first one without any if or ...)




# colocar id no dicionario best ( procurar por abenda valule )






















