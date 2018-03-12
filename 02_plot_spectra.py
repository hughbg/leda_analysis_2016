#!/usr/bin/env python
"""
# 02_plot_spectra.py

Plots the time-averaged spectra.
"""
import matplotlib as mpl
import seaborn as sns
import tables as tb
from leda_cal.leda_cal import *
from leda_cal.skymodel import *

sns.set_style('ticks')
sns.set_context("paper",font_scale=1.5)

def quicklook(filename):
    h5 = tb.open_file(filename)

    T_ant = apply_calibration(h5)
    f_leda = T_ant['f']
    
    print T_ant.keys()
    
    ant_ids = ['252A', '254A', '255A', '254B', '255B']
    pol_id  = 'y'
      
    print("Plotting...")
    fig, ax = plt.subplots(figsize=(8, 6))

    mid = T_ant["252A"].shape[0]/2    # Currently plots the middle bit of the observation
    print T_ant['lst'][mid]           # Print the LST about which we are averaging
    sl  = 50                          # The number of integrations to average either side
    
    print ant_ids[0], mid, sl
    print ant_ids[1], mid, sl
    print ant_ids[2], mid, sl
    
    plt.subplot(2,1,1)
    plt.plot(f_leda, np.median(T_ant[ant_ids[0]][mid-sl:mid+sl], axis=0), label=ant_ids[0])
    plt.plot(f_leda, np.median(T_ant[ant_ids[1]][mid-sl:mid+sl], axis=0), label=ant_ids[1])
    plt.plot(f_leda, np.median(T_ant[ant_ids[2]][mid-sl:mid+sl], axis=0), label=ant_ids[2])
    plt.legend(frameon=False)
    plt.ylabel("Temperature [K]")
    plt.xlim(40, 85)
    plt.minorticks_on()
    plt.ylim(500, 7000)
    
    plt.subplot(2,1,2)
    plt.plot(0, 0)
    plt.plot(f_leda, np.median(T_ant[ant_ids[3]][mid-sl:mid+sl], axis=0), label=ant_ids[3])
    plt.plot(f_leda, np.median(T_ant[ant_ids[4]][mid-sl:mid+sl], axis=0), label=ant_ids[4])
    plt.legend(frameon=False)
    
    plt.xlim(40, 85)
    plt.minorticks_on()
    plt.ylim(500, 7000)
    
    
    plt.xlabel("Frequency [MHz]")
    plt.ylabel("Temperature [K]")
    
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig("figures/compare-spectra.pdf")
    plt.show()
    

if __name__ == "__main__":
    
    import sys
    try:
        filename = sys.argv[1]
    except:
        print "USAGE: ./quicklook.py filename_of_hdf5_observation"
        exit()
    
    quicklook(filename)