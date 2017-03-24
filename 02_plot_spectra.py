#!/usr/bin/env python
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
    
    ant_ids = ['252A', '254A', '255A']
    pol_id  = 'y'
      
    print("Plotting...")
    fig, ax = plt.subplots(figsize=(8, 6))

    mid = T_ant["252A"].shape[0]/2
    sl  = 250
    
    plt.plot(f_leda, np.median(T_ant[ant_ids[0]][mid-sl:mid+sl], axis=0), label=ant_ids[0])
    plt.plot(f_leda, np.median(T_ant[ant_ids[1]][mid-sl:mid+sl], axis=0), label=ant_ids[1])
    plt.plot(f_leda, np.median(T_ant[ant_ids[2]][mid-sl:mid+sl], axis=0), label=ant_ids[2])
    
    plt.xlim(40, 85)
    plt.minorticks_on()
    plt.ylim(500, 7000)
    
    ax.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    
    plt.xlabel("Frequency [MHz]")
    plt.ylabel("Temperature [K]")
    
    plt.legend(frameon=False)
    plt.show()
    

if __name__ == "__main__":
    
    import sys
    try:
        filename = sys.argv[1]
    except:
        print "USAGE: ./quicklook.py filename_of_hdf5_observation"
        exit()
    
    quicklook(filename)