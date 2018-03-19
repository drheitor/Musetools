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

#from multiprocessing import Pool
#from multiprocessing import Process,Queues

#inputs needed to run -atmosphericParameters.dat , radialVelocities.dat names_miles91.lis
#output avg_parameters.txt and banned_stars.txt

FILENAME_OUTPUT_PARAMETERS = 'avg_parameters_coelho.txt'

        
best = {}        


    #open atmosphericParameters and build a best dictionary     
with open('atmosphericParameters.dat') as apar:
   lines = [raw.strip().split() for raw in apar]
   for n, cols in enumerate(lines):
      if len(cols) == 1:
          best[cols[0]] = {'LIB': '', 'id': '', 'rv': ''}
          best[cols[0]]['LIB'] = lines[n+1][0]
          
          #compare em input a id in best
for star in best.keys():
         best[star]['id'] = best[star]['LIB'].strip('star')


with open('radialVelocities.dat') as rvf:
    for cols in [raw.strip().split() for raw in rvf]:
        fname = cols[0].replace('.0', '').replace('-1000', '')
        #check if stars are in radial velocities and vici versa 
        if fname in best.keys():
            best[fname]['rv'] = float(cols[-3])


obj = 'rvTer9'
starpars = {}
with open('atmosphericParameters.dat') as ap:
    for line in [raw.strip().split() for raw in ap]:
        if len(line) == 1:
            star = line[0].split('_')[-1].strip('cb')
            if star not in starpars.keys():
                starpars[star] = {'tpn': [], 'sim': [], 'teff': [], 'lg': [], 'met': [], 'alp': []}
        elif len(line) > 4:
            starpars[star]['tpn'].append(line[0])
            starpars[star]['sim'].append(float(line[1]))
            starpars[star]['teff'].append(float(line[2]))
            starpars[star]['lg'].append(float(line[3]))
            starpars[star]['met'].append(float(line[4]))
            starpars[star]['alp'].append(float(line[5]))


finalpars = {}
rvl = []
mts = []
pfil = open(FILENAME_OUTPUT_PARAMETERS, 'w')
pfil.write('%-8s %5s %8s  %9s   %10s   %10s'%('star_id', 'RV', 'Teff', 'logG', '[Fe/H]', '[Mg/H]\n'))






for s in starpars.keys():
    #print s
    for entry in starpars[s].keys():
        starpars[s][entry] = np.array(starpars[s][entry])
    maskidx = []
    for n, m in enumerate(starpars[s]['tpn']):
        #exclusion criteria
        if (n < 10) and (starpars[s]['lg'][n] < -900.) or (starpars[s]['alp'][n] < -900.) or not ((starpars[s]['lg'][n] <= 4.9 - 0.0002 * starpars[s]['teff'][n]) and (starpars[s]['lg'][n] >= -12.3 + 0.0025 * starpars[s]['teff'][n])):
            maskidx.append(n)
    for entry in starpars[s].keys():
        starpars[s][entry] = np.delete(starpars[s][entry], maskidx)
    if len(starpars[s]['teff']) > 1:
        finalpars[s] = {}
        weight = 1/(np.array(starpars[s]['sim'])**2)
        teff = np.average(starpars[s]['teff'], weights=weight)
        finalpars[s]['teff'] = teff
        logg = np.average(starpars[s]['lg'], weights=weight)
        finalpars[s]['logg'] = logg
        met = np.average(starpars[s]['met'], weights=weight)
        finalpars[s]['met'] = met
        mts.append(met)
        alp =np.average(starpars[s]['alp'], weights=weight)
        finalpars[s]['alp'] = alp
        eteff = np.sqrt(np.average((starpars[s]['teff']-teff)**2, weights=weight))
        finalpars[s]['eteff'] = eteff
        elogg = np.sqrt(np.average((starpars[s]['lg']-logg)**2, weights=weight))
        finalpars[s]['elogg'] = elogg
        emet = np.sqrt(np.average((starpars[s]['met']-met)**2, weights=weight))
        finalpars[s]['emet'] = emet
        ealp = np.sqrt(np.average((starpars[s]['alp']-alp)**2, weights=weight))
        finalpars[s]['ealp'] = ealp
        spname = s+'cb'
        rv = best[spname]['rv']
        rvl.append(rv)
        finalpars[s]['rv'] = rv
        pfil.write('%-8s %7.2f  %4.0f±%-4.0f %5.2f±%4.2f %6.2f±%4.2f %6.2f±%4.2f\n'%(s, rv, teff, eteff, logg, elogg, met, emet, alp, ealp))
    else:
        print("**Attention** Star '{}' will not be saved in parameters file".format(s))
        ban = open("banned_stars.txt","a") 
        ban.write("**Attention** Star '{}' will not be saved in parameters file \n".format(s) ) 
        ban.close() 
pfil.close()

print(" Saved in file '{}'".format(FILENAME_OUTPUT_PARAMETERS))

