#!/usr/bin/env python
"""
# 11_plot_nw.py

Plots the noise wave as computed from calibration

"""
import matplotlib as mpl
import seaborn as sns
import tables as tb
from leda_cal.leda_cal import *
from leda_cal.skymodel import *
from leda_cal.dpflgr import *

from lmfit import minimize, Parameters, report_fit

sns.set_style('ticks')
sns.set_context("paper",font_scale=1.5)

def quicklook(filename):
    
    balun_loss = hkl.load('cal_data/balun_loss.hkl')
    vna_cal    = hkl.load('cal_data/vna_calibration.hkl')
    
    ant_id = 'a252x'
    ra = vna_cal[ant_id]["ra"]
    rl = vna_cal[ant_id]["rl"]
    f  = vna_cal["f"]
    F0 = compute_F(ra, rl)
    G = compute_G(rl)
    
    T_noise = compute_noisewave('a252x')
    
    plt.plot(f, T_noise)
    plt.show()
    

if __name__ == "__main__":
    
    import sys
    try:
        filename = sys.argv[1]
    except:
        print "USAGE: ./quicklook.py filename_of_hdf5_observation"
        exit()
    
    quicklook(filename)