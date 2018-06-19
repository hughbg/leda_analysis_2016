#!/usr/bin/env python
"""
# 04_compare_spectra.py

Compare spectra across antennas by dividing through by average
to see fractional differences

"""

import seaborn as sns
import tables as tb
from leda_cal.skymodel import *
from leda_cal.leda_cal import *
from leda_cal.git import get_repo_fingerprint
from leda_cal.dpflgr import simple_flag
sns.set_style('ticks')
sns.set_context("paper",font_scale=1.7)

def quicklook(filename):
    h5 = tb.open_file(filename)

    T_ant = apply_calibration(h5)
    f_leda = T_ant['f']
    
    ant_ids = ['252A', '254A', '255A', '254B', '255B']
    pol_id  = 'y'
      
    print("Plotting...")
    fig, axes = plt.subplots(figsize=(8, 6))
    #plt.suptitle(h5.filename)
    
    lst_stamps = T_ant['lst']
    utc_stamps = T_ant['utc']
    xlims = (f_leda[0], f_leda[-1])
    #ylims = mdates.date2num((T_ant['utc'][0], T_ant['utc'][-1]))
    #hfmt = mdates.DateFormatter('%m/%d %H:%M')
    ylims = (T_ant['lst'][0], T_ant['lst'][-1])
    
    T_avg = (T_ant["252A"] + T_ant["254A"] + T_ant["255A"]) / 3
    
    mid = T_ant["252A"].shape[0]/2
    sl  = 250
    
    plt.subplot(2,1,1)
    T_avg = np.median(T_avg[mid-sl:mid+sl], axis=0)
    T_252 = np.median(T_ant[ant_ids[0]][mid-sl:mid+sl], axis=0) - T_avg
    T_254 = np.median(T_ant[ant_ids[1]][mid-sl:mid+sl], axis=0) - T_avg
    T_255 = np.median(T_ant[ant_ids[2]][mid-sl:mid+sl], axis=0) - T_avg
    
    T_252 = simple_flag(T_252, 20.0)
    T_254 = simple_flag(T_254, 20.0)
    T_255 = simple_flag(T_255, 20.0)
    
    T_252 = T_252 / T_avg
    T_254 = T_254 / T_avg
    T_255 = T_255 / T_avg
    
    plt.plot(f_leda, T_252, label='252A')
    plt.plot(f_leda, T_254, label='254A')
    plt.plot(f_leda, T_255, label='255A')
    
    plt.xlabel("Frequency [MHz]")
    plt.ylabel("Fractional difference")

    plt.xlim(40, 85)
    plt.ylim(-0.07, 0.07)
    plt.legend(loc=4, frameon=True, ncol=3)
    plt.minorticks_on()
    
    plt.subplot(2,1,2)
    T_avg = (T_ant["254B"] + T_ant["255B"]) / 2.0
    
    mid = T_ant["252A"].shape[0]/2
    sl  = 250
    
    T_avg = np.median(T_avg[mid-sl:mid+sl], axis=0)
    #T_252 = np.median(T_ant[ant_ids[0]][mid-sl:mid+sl], axis=0) - T_avg
    T_254 = np.median(T_ant[ant_ids[3]][mid-sl:mid+sl], axis=0) - T_avg
    T_255 = np.median(T_ant[ant_ids[4]][mid-sl:mid+sl], axis=0) - T_avg
    
    #T_252 = simple_flag(T_252, 20.0)
    T_254 = simple_flag(T_254, 20.0)
    T_255 = simple_flag(T_255, 20.0)
    
    #T_252 = T_252 / T_avg
    T_254 = T_254 / T_avg
    T_255 = T_255 / T_avg
    
    #plt.plot(f_leda, T_252, label='252A')
    plt.plot(0,0)
    plt.plot(f_leda, T_254, label='254B')
    plt.plot(f_leda, T_255, label='255B')
    
    plt.xlabel("Frequency [MHz]")
    plt.ylabel("Fractional difference")

    plt.xlim(40, 85)
    plt.ylim(-0.07, 0.07)
    plt.legend(loc=4, frameon=True, ncol=3)
    plt.minorticks_on()
    
    plt.text(0.005, 0.005, get_repo_fingerprint(), transform=fig.transFigure, size=8)
    
    plt.savefig("figures/compare-ants.pdf")
    plt.show()
    
    

if __name__ == "__main__":
    
    import sys
    try:
        filename = sys.argv[1]
    except:
        print "USAGE: ./quicklook.py filename_of_hdf5_observation"
        exit()
    
    quicklook(filename)
