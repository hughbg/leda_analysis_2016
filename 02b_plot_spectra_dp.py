#!/usr/bin/env python
"""
# 02b_plot_spectra_dp.py

Plot spectra for single antenna in dual panels:
1) Spectra at LST where power is minimal
2) Power vs time at 60 MHz

"""
import matplotlib as mpl
import seaborn as sns
import tables as tb
from leda_cal.leda_cal import *
from leda_cal.skymodel import *

from leda_cal.dpflgr import *

sns.set_style('ticks')
sns.set_context("paper",font_scale=1.5)

def quicklook(filename):
    h5 = tb.open_file(filename)

    T_ant = apply_calibration(h5)
    f_leda = T_ant['f']
    lst    = T_ant['lst']
    
    print T_ant.keys()
    
    ant_ids = ['252A']
      
    print("Plotting...")
    fig, ax = plt.subplots(figsize=(8, 6))

    mid = closest(lst, 11)
    print T_ant['lst'][mid]
    sl  = 20
    
    T_flagged = rfi_flag(T_ant[ant_ids[0]], freqs=f_leda)
    
    plt.title("252A")
    plt.subplot(2,1,1)
    plt.plot(f_leda, T_flagged[mid-sl:mid+sl].mean(axis=0), label='LST 11:00', c='#333333')
    plt.axvline(x=60, ls='dashed', c='#888888')
    plt.xlim(40, 85)
    plt.minorticks_on()
    plt.ylim(500, 7000)
    
    ax.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    
    plt.xlabel("Frequency [MHz]")
    plt.ylabel("Temperature [K]")
    plt.legend(frameon=False)
    
    plt.subplot(2,1,2)
    sl = 10
    #mid = closest(f_leda, 40)
    #plt.plot(lst, T_flagged[:, mid-sl:mid+sl].mean(axis=1), label='40 MHz')
    mid = closest(f_leda, 60)
    plt.plot(lst, T_flagged[:, mid-sl:mid+sl].mean(axis=1), label='60 MHz', c='#333333')
    plt.axvline(x=11, ls='dashed', c='#888888')
    #mid = closest(f_leda, 80)
    #plt.plot(lst, T_flagged[:, mid-sl:mid+sl].mean(axis=1), label='80 MHz')

    #plt.xlim(40, 85)
    plt.xlabel("Sidereal time [hr]")
    plt.ylabel("Temperature [K]")
    plt.minorticks_on()
    plt.xlim(0, 24)
    plt.xticks([0, 4, 8, 12, 16, 20])
    
    ax.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    
    
    plt.legend(frameon=False, loc=2)
    plt.tight_layout()
    plt.savefig("figures/a252-cuts.pdf")
    plt.show()
    sl  = 250
    plt.plot(f_leda, T_flagged[mid-sl:mid+sl].mean(axis=0) - T_flagged.mean(axis=0))
    plt.show()
    

if __name__ == "__main__":
    
    import sys
    try:
        filename = sys.argv[1]
    except:
        print "USAGE: ./quicklook.py filename_of_hdf5_observation"
        exit()
    
    quicklook(filename)
