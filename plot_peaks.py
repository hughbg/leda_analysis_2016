#!/usr/bin/env python
"""
# 02_plot_spectra.py

Plots the time-averaged spectra.
"""
import matplotlib as mpl
import os
import seaborn as sns
import tables as tb
from leda_cal.leda_cal import *
from leda_cal.skymodel import *
from leda_cal.dpflgr import *
from leda_cal.git import get_repo_fingerprint

sns.set_style('ticks')
sns.set_context("paper",font_scale=1.5)

def quicklook(filename, save, flag, noshow):
    h5 = tb.open_file(filename)

    T_ant = apply_calibration(h5)
    f_leda = T_ant['f']
    
    ant_ids = ['252A', '254A', '255A', '252B', '254B', '255B'] # Added 252B, why was it missing? HG
    mmax = 0
    for ant in ant_ids:
      if flag:
        T_flagged = rfi_flag(T_ant[ant], freqs=f_leda)
        T_ant[ant] = np.ma.filled(T_flagged, 0)

      if flag: tmax = np.max(T_ant[ant])
      else: tmax = np.percentile(T_ant[ant], 99.95)
      if tmax  > mmax: mmax = tmax  

    print("Plotting...")
    fig, ax = plt.subplots(figsize=(8, 6))

    mid = T_ant["252A"].shape[0]/2    # Currently plots the middle bit of the observation
    print T_ant['lst'][mid]           # Print the LST about which we are averaging
    sl  = 50                          # The number of integrations to average either side
    
    print ant_ids[0], mid, sl
    print ant_ids[1], mid, sl
    print ant_ids[2], mid, sl
    
    plt.subplot(2,1,1)
    plt.plot(f_leda, np.max(T_ant[ant_ids[0]], axis=0), label=ant_ids[0])
    plt.plot(f_leda, np.max(T_ant[ant_ids[1]], axis=0), label=ant_ids[1])
    plt.plot(f_leda, np.max(T_ant[ant_ids[2]], axis=0), label=ant_ids[2])
    plt.legend(frameon=False)
    plt.ylabel("Temperature [K]")
    plt.xlim(30, 85)
    plt.minorticks_on()
    plt.ylim(500, mmax)
    
    plt.subplot(2,1,2)
    plt.plot(0, 0)
    plt.plot(f_leda, np.max(T_ant[ant_ids[3]], axis=0), label=ant_ids[3])
    plt.plot(f_leda, np.max(T_ant[ant_ids[4]], axis=0), label=ant_ids[4])
    plt.plot(f_leda, np.max(T_ant[ant_ids[4]], axis=0), label=ant_ids[5])

    plt.legend(frameon=False)
    
    plt.xlim(30, 85)
    plt.minorticks_on()
    plt.ylim(500, mmax)
    
    
    plt.xlabel("Frequency [MHz]")
    plt.ylabel("Temperature [K]")
    
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.text(0.005, 0.005, get_repo_fingerprint(), transform=fig.transFigure, size=8)
    plt.savefig("figures/compare-spectra.pdf")
    
    if save:
        if flag:
            plt.savefig("peaks_"+os.path.basename(filename)[:-3]+"_"+"flagged.png")
        else:
            plt.savefig("peaks_"+os.path.basename(filename)[:-3]+".png")
    if not noshow:
        plt.show()

if __name__ == "__main__":
    import optparse, sys

    usage = '%prog [opts] filename_of_hdf5_observation'
    o = optparse.OptionParser()
    o.set_usage(usage)
    o.set_description(__doc__)
    o.add_option('--flag', dest='flag', action='store_true', default=False,
                 help='Apply flagging. Default: False')
    o.add_option('--save', dest='save', action='store_true', default=False,
                 help="Save the plot to an image file, with filename the same as the h5 but png extension. Default: False.")
    o.add_option('--noshow', dest='noshow', action='store_true', default=False,
                 help="Don't display the plot on screen. Useful for batch runs. Default: False.")

    opts, args = o.parse_args(sys.argv[1:])


    if len(args) != 1:
        o.print_help()
        exit(1)
    else:
        filename = args[0]
        
    quicklook(filename, opts.save, opts.flag, opts.noshow)
