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

# MILES=5, COELHO=4
LIB = 4

FINAL_OUTPUT = 'Avg_par.dat'

if (LIB==4):
        FINAL_OUTPUT = 'Avg_par_Coelho.dat'

if (LIB==5):
        FINAL_OUTPUT = 'Avg_par_Miles.dat'        

################################################## MAIN #####################################################

#defines lim 
lwlim = 4200
uwlim = 6000


################################################## FLAG.DAT ###########################################################

# this function builds the flag.dat file containing etoile input parameters
def genflag(type):
    
    #this function builds tha flag.dat file containing etoile parameters
    global lwlim
    global uwlim
    global LIB
    parset = {'RV': 0, 'PAR': 0, 'COMP': 0}
    parset[type] = 1
    if LIB == 5:
        libpath = '../code/dnm/parameters_standard_V20130201.dat'
    elif LIB == 4:
        libpath = '../code/dnm/parameters_standard_coelho05.dat'
    with open('flag.dat', 'w') as flgfil:
        flgfil.write('DERIVE_RADIAL_VELOCITY....... '+str(parset['derrv'])+
                     '\n     write_logs................... 2\n \nNormalise_Object_Continuum... 1'
                     '\n     display_spectrum............. 0\n     reject_high.................. 0'
                     '\n     degree....................... 4\n     sigma_rej.................... 1.8'
                     '\n     nb_kept_pixel_min............ 15\n \nNormalise_Object_Fluxes...... 1\n '
                     '\nSelect_Wavelength_Domain_Temp 1\n     wave_blue.................... '+'%.2f'%lwlim+'\n     '
                     'wave_red..................... '+'%.2f'%uwlim+'\n \nNormalise_Template_Continuum. 1\n     '
                     'display_spectrum............. 0\n     reject_high.................. 0\n     '
                     'degree....................... 4\n     sigma_rej.................... 1.8\n     '
                     'nb_kept_pixel_min............ 15\n \nNormalise_Template_Fluxes.... 1\n \n'
                     'Prepare_Both_Spectra......... 0\n     vmax......................... 700.0\n \n'
                     'Calculate_Analyse_CC_Peak\n     nbst0........................  60\n     '
                     'nbst1........................ 225\n     step0........................  10.0\n     '
                     'step1........................   0.2\n \nDisplay_Overplotted_Spectra.. 1\n \n     '
                     '-----------------------------------\n \nLIBRARY_OF_ATMOSPHERIC_PARAMETERS\n'
                     'Param_Std '+libpath+'\n\n     '
                     '-----------------------------------\n \nCOMPARE_SPECTRA.............. '+
                     str(parset['compare'])+'\nDERIVE_ATMOSPHERIC_PARAMETERS '+
                     str(parset['deratmos'])+'\n     Match_Orders................. 0\n     '
                     'Convolve..................... 0\n     Graph........................ 0\n     '
                     'Write........................ 1\n     Number_Of_Degrees_Of_Freedom. 4')

####################################################################
        
 def pick_best():
     
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
         
    return best

################################################## RV ###########################################################



###################################################ATMOSPHERICPARAMETERS###################################################


###################################################CHECK FIT##############################################################

#checking etoile fit

for _spname in inlist:
    spname = _spname.split('.')[0]
    sp_c = np.genfromtxt(spname+'_fit.tab', unpack=True, skip_header=1)
    plt.figure(figsize=(16,8))
    plt.plot(sp_c[0], sp_c[1], 'k-', lw=0.7, label='observed')
    plt.plot(sp_c[0], sp_c[2], 'r-', lw=0.7, label='library')
    plt.gca().set_xlabel('Wavelength [$\AA$]')
    plt.gca().set_ylabel('normalized flux')
    plt.axis([5000, 5900, 0.0, 2.2])
    plt.legend()
    plt.tight_layout()
    plt.savefig("compared-{}.pdf".format(spname))


################################################### GET OUTPUT ##############################################################


obj = 'Ter9'
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
for s in sorted(starpars.keys()):
    #print s
    for entry in starpars[s].keys():
        starpars[s][entry] = np.array(starpars[s][entry])
    maskidx = []
    starpars[s]['sim'] = starpars[s]['sim']/starpars[s]['sim'][0]
    for n, m in enumerate(starpars[s]['sim']):
        #exclusion criteria
        if (m > 1.1) or (starpars[s]['lg'][n] < -900.) or (starpars[s]['alp'][n] < -900.) or not (        (starpars[s]['lg'][n] <= 4.9 - 0.0002 * starpars[s]['teff'][n]) and        (starpars[s]['lg'][n] >= -12.3 + 0.0025 * starpars[s]['teff'][n])):
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
        spname = obj+'_'+s+'cb'
        rv = bestmatch[spname]['rv']
        rvl.append(rv)
        finalpars[s]['rv'] = rv
        pfil.write('%-8s %7.2f  %4.0f±%-4.0f %5.2f±%4.2f %6.2f±%4.2f %6.2f±%4.2f\n'%(s, rv, teff, eteff, logg, elogg, met, emet, alp, ealp))
    else:
        print("**Attention** Star '{}' will not be saved in parameters file".format(s))

pfil.close()

print("Just saved file '{}'".format(FILENAME_OUTPUT_PARAMETERS))

