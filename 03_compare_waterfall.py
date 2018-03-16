#!/usr/bin/env python
"""
# 03_compare_waterfall.py

Compare waterfall plots across antennas (subtract average & plot)

"""
import seaborn as sns
import tables as tb
from leda_cal.skymodel import *
from leda_cal.leda_cal import *
sns.set_style('white')
sns.set_context("poster",font_scale=.75)

def quicklook(filename):
    h5 = tb.open_file(filename)

    T_ant = apply_calibration(h5)
    f_leda = T_ant['f']
    
    ant_ids = ['252A', '254A', '255A']
    pol_id  = 'y'
      
    print("Plotting...")
    fig, axes = plt.subplots(figsize=(12, 6))
    #plt.suptitle(h5.filename)
    
    lst_stamps = T_ant['lst']
    utc_stamps = T_ant['utc']
    xlims = (f_leda[0], f_leda[-1])
    #ylims = mdates.date2num((T_ant['utc'][0], T_ant['utc'][-1]))
    #hfmt = mdates.DateFormatter('%m/%d %H:%M')
    ylims = (T_ant['lst'][0], T_ant['lst'][-1])
    
    T_avg = (T_ant["252A"] + T_ant["254A"] + T_ant["255A"]) / 3


    for ii, key in enumerate(ant_ids):
        ax = fig.add_subplot(1, 3, ii+1)
        
        T_flagged = T_ant[key]
        #T_flagged = rfi_flag(T_ant[key], thr_f=0.07, thr_t=0.07, rho=1.5,
        #         bp_window_f=8, bp_window_t=8, 
        #         max_frac_f=0.3, max_frac_t=0.3)
        
        T_ant[key] = T_flagged
        
        im = plt.imshow(T_flagged - T_avg, # / np.median(xx, axis=0), 
                   cmap='viridis', aspect='auto',
                   interpolation='nearest',
                   clim=(-250, 250),
                   extent=(xlims[0], xlims[1], ylims[1], ylims[0])
                   )
        plt.title(ant_ids[ii])
        plt.xlabel("Frequency [MHz]")
        #ax.yaxis_date()
        #ax.yaxis.set_major_formatter(hfmt)
        #
    
    plt.subplot(1,3,1)
    plt.ylabel("LST [hr]")
    
    fig.subplots_adjust(left=0.06)
    fig.subplots_adjust(right=0.875)
    cbar_ax = fig.add_axes([0.9, 0.125, 0.025, 0.75])
    cbar = fig.colorbar(im, cax=cbar_ax)
    
    plt.subplot(1,3,3)
    #cbar = plt.colorbar()
    cbar.set_label("Temperature [K]")
    
    #plt.tight_layout()
    plt.show()
    exit()	# I get an error if it goes further
    plt.xlabel("Frequency [MHz]")

    plt.ylabel("LST [hr]")
    
    #fig.subplots_adjust(left=0.06)
    #fig.subplots_adjust(right=0.875)
    #cbar_ax = fig.add_axes([0.9, 0.125, 0.025, 0.75])
    #cbar = fig.colorbar(im, cax=cbar_ax)
    cbar = plt.colorbar()
    
    #plt.subplot(1,3,3)
    #cbar = plt.colorbar()
    cbar.set_label("Temperature [K]")
    
    #plt.tight_layout()
    plt.show()
    
    

if __name__ == "__main__":
    
    import sys
    try:
        filename = sys.argv[1]
    except:
        print "USAGE: ./quicklook.py filename_of_hdf5_observation"
        exit()
    
    quicklook(filename)
