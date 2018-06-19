#!/usr/bin/env python
"""
# 06_plot_rfi.py

Plots the RFI flags

"""
import numpy as np
import seaborn as sns
import tables as tb
from leda_cal.skymodel import *
from leda_cal.leda_cal import *
from leda_cal.dpflgr import *

from scipy.stats import kurtosis, scoreatpercentile as percentile

sns.set_style('ticks')
sns.set_context("paper",font_scale=1.5)

def quicklook(filename, flatten, ant='252A'):
    h5 = tb.open_file(filename)

    T_ant = apply_calibration(h5)
    f_leda = T_ant['f']
    
    ant_ids = [ant,]
      
    print("Plotting %s..." % ant_ids[0])
    fig, axes = plt.subplots(figsize=(12, 6), nrows=1, ncols=1)
    #plt.suptitle(h5.filename)
    
    lst_stamps = T_ant['lst']
    utc_stamps = T_ant['utc']
    xlims = (f_leda[0], f_leda[-1])
    #ylims = mdates.date2num((T_ant['utc'][0], T_ant['utc'][-1]))
    #hfmt = mdates.DateFormatter('%m/%d %H:%M')
    ylims = (T_ant['lst'][0], T_ant['lst'][-1])
    T_flagged = T_ant[ant_ids[0]]
    #T_flagged = np.fft.fft(T_flagged, axis=0)
    #T_flagged -= T_flagged.mean(axis=0)
    #T_flagged = 10*np.log10(np.abs(np.fft.ifft(T_flagged)))

    T_flagged = rfi_flag(T_flagged, freqs=f_leda)
    
    if flatten:
        abp = np.ma.median(T_flagged.data, axis=0)
        abp /= np.ma.median(abp)
        T_flagged /= abp
    clim = (percentile(T_flagged.compressed(), 5), percentile(T_flagged.compressed(), 95))
        
    im = plt.imshow(T_flagged, # / np.median(xx, axis=0), 
                    cmap='jet', aspect='auto',
                    interpolation='nearest',
                    clim=clim,
                    extent=(xlims[0], xlims[1], ylims[1], ylims[0]))
    plt.title(ant_ids[0])
    plt.xlabel("Frequency [MHz]")

    plt.ylabel("LST [hr]")
    plt.colorbar()
    plt.savefig("figures/rfi-flagged.pdf")
    plt.show()
    
    plt.figure()
    #plt.plot(f_leda, np.sum(T_flagged.mask, axis=0).astype('float') / T_flagged.mask.shape[0], label='total')
    day = T_flagged[0:2000].mask
    night = T_flagged[2250:2750].mask
    plt.plot(f_leda, np.sum(night, axis=0).astype('float') / night.shape[0], label='night')
    plt.plot([0])
    plt.plot(f_leda, np.sum(day, axis=0).astype('float') / day.shape[0], label='day')
    plt.xlim(40, 85)
    plt.ylim(-0.025, 0.25)
    
    plt.title(ant_ids[0])
    plt.xlabel("Frequency [MHz]")
    plt.ylabel("Flagged fraction")
    plt.minorticks_on()
    plt.legend(frameon=True, loc=2)
    plt.tight_layout()
    plt.savefig("figures/rfi-fraction.pdf")
    plt.show()
    
    plt.figure()
    plt.plot(f_leda, kurtosis(T_flagged, axis=0))
    plt.title(ant_ids[0])
    plt.ylabel("Kurtosis")
    plt.xlabel("Frequency [MHz]")
    plt.xlim(40, 85)
    plt.ylim(-50, 1600)
    plt.minorticks_on()
    plt.show()
    
    plt.figure()
    

if __name__ == "__main__":
    import optparse, sys

    usage = '%prog [opts] filename_of_hdf5_observation'
    o = optparse.OptionParser()
    o.set_usage(usage)
    o.set_description(__doc__)
    o.add_option('--ant', dest='ant', action='store', default='252A', 
                 help='Name of the antenna to plot. Default: 252A')
    o.add_option('--flatten', dest='flatten', action='store_true', default=False,
                 help='Apply a crude bandpass derived from the data. Default: False')
    
    opts, args = o.parse_args(sys.argv[1:])
    
    if len(args) != 1:
        o.print_help()
        exit(1)
    else:
        filename = args[0]
        
    quicklook(filename, opts.flatten, ant=opts.ant)
