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
#import multiprocessing as mp
import codecs
import sys


DEFAULT_MIDDLE = "coelho"

FILENAMES_OUTPUT = ["plot-alpha.pdf", "plot-hr.pdf", "plot-distribution.pdf"]




def read_parameters_file(filename):

    """
    
    Returns:
        dict: {"star_id": (dictionary with keys 'rv', 'teff', 'eteff', 'logg', 'elogg', 'met', 'emet', 'alp', 'ealp'), ...}   

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
        pfil.write('%-8s %7.2f  %4.0f±%-4.0f %5.2f±%4.2f %6.2f±%4.2f %6.2f±%4.2f\n'%(s, rv, teff, eteff, logg, elogg, met, emet, alp, 



star_id     RV     Teff       logG       [Fe/H]      [Mg/H]
0056       10.02  4844±344   2.17±0.87  -1.26±0.62   0.40±0.15
0072       62.72  4591±187   2.24±0.52  -0.63±0.28   0.39±0.20
0077      -78.41  4413±126   2.13±0.39   0.06±0.15  -0.01±0.16
0081       44.42  4573±188   2.17±0.56  -0.58±0.35   0.34±0.26
0084       38.60  4557±175   2.06±0.54  -0.73±0.47   0.38±0.20
0086      424.29  5028±453   2.20±0.86  -1.76±0.48   0.41±0.17
0089       61.45  4616±267   2.17±0.68  -0.53±0.56   0.22±0.23
0         1         2         3         4         5         6
0123456789012345678901234567890123456789012345678901234567890123456789
    """

    finalpars = {}
    with codecs.open(filename,'r',encoding='utf8') as file:
        
        for i, line in enumerate(file):
            if i == 0:
                # Skips first line
                continue
        
            s = str(line[0:4])
        
            row = finalpars[s] = {}
            
#            print("----------->{}".format(line))
#            print("----------->0123456789012345678901234567890123456789012345678901234567890123456789           ")
            row["rv"] = float(line[5:16])
            row["teff"] = float(line[18:22])
            row["eteff"] = float(line[23:26])
            row["logg"] = float(line[27:33])
            row["elogg"] = float(line[34:38])
            row["met"] = float(line[39:45])
            row["emet"] = float(line[46:50])
            row["alp"] = float(line[51:57])
            row["ealp"] = float(line[58:62])

    return finalpars


def plot_plots(finalpars):

    rvmin = -100
    rvmax = 100

    rvl = [row["rv"] for row in finalpars.values()]
    mts = [row["met"] for row in finalpars.values()]
   
    rvl = np.array(rvl)
    delmask = []
    for n in range(len(rvl)):
        if (rvl[n] < rvmin) or (rvl[n] > rvmax):
            delmask.append(n)
    rvl = np.delete(rvl, delmask)
    mnrv = abs(min(rvl))
    nrvl = rvl + mnrv
    mxrv = max(nrvl)



    f1 = plt.figure(figsize = (10,6))
    grid = gs.GridSpec(1, 2, width_ratios=[10, 0.5]) 
    ax1 = plt.subplot(grid[0])
    cmap = plt.get_cmap('plasma')
    for s in sorted(finalpars.keys()):
        if (finalpars[s]['rv'] > rvmin) and (finalpars[s]['rv'] < rvmax):
            ax1.plot(finalpars[s]['met'], finalpars[s]['alp'], marker='.', color=cmap((finalpars[s]['rv']+mnrv)/mxrv))
    ax1.set_xlabel('[Fe/H]')
    ax1.set_ylabel('[Mg/Fe]')
    cbr = plt.subplot(grid[1])
    cmap = mpl.cm.plasma
    norm = mpl.colors.Normalize(vmin=min(rvl), vmax=max(rvl))
    cb1 = mpl.colorbar.ColorbarBase(cbr, cmap=cmap,
                                    norm=norm,
                                    orientation='vertical')
    cb1.set_label('Radial Velocity')
    plt.tight_layout()
    plt.savefig(FILENAMES_OUTPUT[0])

    # =========================
    f2 = plt.figure(figsize=(7,6))
    grid2 = gs.GridSpec(1, 2, width_ratios=[10, 0.5]) 
    ax2 = plt.subplot(grid2[0])
    for s in sorted(finalpars.keys()):
        if (finalpars[s]['rv'] > -100) and (finalpars[s]['rv'] < 100):
            ax2.plot(finalpars[s]['teff'], finalpars[s]['logg'], marker='.', color=cmap((finalpars[s]['rv']+mnrv)/mxrv))
    ax2.invert_xaxis()
    ax2.invert_yaxis()
    ax2.set_xlabel('Teff')
    ax2.set_ylabel('logG')
    cbr2 = plt.subplot(grid2[1])
    cb2 = mpl.colorbar.ColorbarBase(cbr2, cmap=cmap,
                                    norm=norm,
                                    orientation='vertical')
    cb2.set_label('Radial Velocity')
    plt.savefig(FILENAMES_OUTPUT[1])

    # =========================
    f3 = plt.figure(figsize=(10,6))
    ax3 = f3.add_subplot(111)
    density = gaussian_kde(mts)
    xs = np.linspace(min(mts), max(mts), 200)
    density.covariance_factor = lambda : .12
    density._compute_covariance()
    ax3.plot(xs, density(xs), 'b-')
    ax3.set_xlabel('[Fe/H]')
    ax3.set_ylabel('Density')
    plt.savefig(FILENAMES_OUTPUT[2])
        


    # In[ ]: __name__ == "__main__":

    if len(sys.argv) > 1:
        middle = sys.argv[1]
    else:
        print("Middle of filename not specified, using default '{}'".format(DEFAULT_MIDDLE))
        middle = DEFAULT_MIDDLE
           
    filename_input = "parameters_{}.dat".format(middle)
    
    print("Reading file '{}'...".format(filename_input))
  

    finalpars = read_parameters_file(filename_input)
    plot_plots(finalpars)
    
    print("*** Output files: ***")
    for filename in FILENAMES_OUTPUT:
        print("    - {}".format(filename))
    

