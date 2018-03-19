#!/usr/bin/env python
# coding: utf-8


# import libraries
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec as gs
import subprocess as sp
from scipy.stats import gaussian_kde

LIB = 4


def get_index():
    
    bestmatch = {}
    with open('atmosphericParameters.dat') as apar:
        lines = [raw.strip().split() for raw in apar]
        for n, cols in enumerate(lines):
            if len(cols) == 1:
                bestmatch[cols[0]] = {'LIB': '', 'id': '', 'rv': ''}
                bestmatch[cols[0]]['LIB'] = lines[n+1][0]
                
        for star in bestmatch.keys():
            bestmatch[star]['id'] = bestmatch[star]['LIB'].strip('star')
    
    return bestmatch


bestmatch = get_index()

with open('./templatesForCC.txt', 'w') as tplf:
    tplf.write(str(len(bestmatch.keys()))+'\n')
    for n, star in enumerate(bestmatch.keys()):
        tplf.write(star+' '+str(n+1)+' '+bestmatch[star]['LIB']+' '+bestmatch[star]['id']+'\n')    















