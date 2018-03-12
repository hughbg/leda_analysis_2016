#!/usr/bin/env python

"""
plot_spa.py

Script for plotting SPA files from anritsu VNA
"""

import os, sys
import numpy as np
import pylab as plt

def read_spa(filename):
    """ Read an SPA file """
    
    ff = open(filename, 'r')
    dd = []
    
    # Read file
    lines = ff.readlines()
    
    # find last comment line (this is row headers)
    data = {}
    d_id = 0
    for line in lines:
        if line.startswith('# Begin TRACE'):
            data[d_id] = []
        if line.startswith('P_'):
            dpoint = map(float, line[2:].replace('=', ' ').replace(',', '').replace('MHz', '').split())
            data[d_id].append(dpoint)
        if line.startswith('# Data Done'):
            d_id += 1
    
    return [np.array(dd) for dd in data.values()]
    
    
def plot_spa(filename, save_fig=True, show_fig=True, c='#cc0000', label=None):
    """ Plot S2P file """
    
    data = read_spa(filename)
    
    ff  = data[0][:, 2]
    tr1 = data[0][:, 1]
    tr2 = data[1][:, 1]
    
    #plt.figure("SPA", figsize=(8,6))
    #plt.title(os.path.split(filename)[-1])
    ## FIX FLIPPED DATA WHERE TR1 <-> TR2
    plt.plot(ff, tr1, c=c, label=label)
    plt.xlabel("Frequency [MHz]")
    plt.ylabel("Power [dB]")
    plt.minorticks_on()
    
    if save_fig:
        plt.savefig(os.path.splitext(filename)[0] + '.pdf')
    if show_fig:
        plt.show()   
    else:
        #plt.clf() 
        pass
  

if __name__ == "__main__":
    try:
        fname = sys.argv[1]
        #fname = '2b2a-rj45-Ncal_allTERMoutput.s2p'
    except:
        print "USAGE: ./plot_spa.py FILENAME"
        exit()
    
    plot_spa(fname)