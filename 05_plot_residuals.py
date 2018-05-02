#!/usr/bin/env python
"""
# 05_plot_residuals.py

Fits a polynomial and subtracts it, then plots result.

"""
import matplotlib as mpl
import seaborn as sns
import tables as tb
from leda_cal.leda_cal import *
from leda_cal.skymodel import *
from leda_cal.dpflgr import *

sns.set_style('ticks')
sns.set_context("paper",font_scale=1.5)

def fit_poly(x, y, n=5, log=True):
    """ Fit a polynomial to x, y data 
    
    x (np.array): x-axis of data (e.g. frequency)
    y (np.array): y-axis of data (e.g temperature)
    n (int): number of terms in polynomial (defaults to 5)
    """
    
    x_g = x 
    x = np.ma.array(x, mask=y.mask).compressed()
    y = y.compressed()
    
    if log:
        print "HERE"
        yl = np.log10(y)
    else:
        yl = y
    fit = np.polyfit(x, yl, n)
    print fit
    p = np.poly1d(fit)
    
    if log:
        return 10**(p(x_g))
    else:
        return p(x_g)


def fit_poly_log(x, y, n=5):
    """ Fit a polynomial to x, y data 
    
    x (np.array): x-axis of data (e.g. frequency)
    y (np.array): y-axis of data (e.g temperature)
    n (int): number of terms in polynomial (defaults to 5)
    """

    x_g = x 
    x = np.ma.array(x, mask=y.mask).compressed()
    y = y.compressed()
    yl = np.log10(y)
    fit = np.polyfit(x, yl, n)
    print fit
    p = np.poly1d(fit)
    return 10**(p(x_g))

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
    
    
    
    d0 = T_ant[ant_ids[0]][mid-sl:mid+sl]
    d1 = T_ant[ant_ids[1]][mid-sl:mid+sl]
    d2 = T_ant[ant_ids[2]][mid-sl:mid+sl]
    
    d0 = rfi_flag(d0, freqs=f_leda)
    d0 = np.median(d0, axis=0)

    d1 = rfi_flag(d1, freqs=f_leda)
    d1 = np.median(d1, axis=0)
    
    d2 = rfi_flag(d2, freqs=f_leda)
    d2 = np.median(d2, axis=0)
    
    #plt.imshow(d0, cmap='viridis', aspect='auto')
    #plt.colorbar()
    #plt.show()
    #exit()
    
    
    
    print d0
    
    plt.plot(f_leda, d0 - fit_poly_log(f_leda, d0, 5), label=ant_ids[0])
    plt.plot(f_leda, d1 - fit_poly_log(f_leda, d1, 5), label=ant_ids[1])
    plt.plot(f_leda, d2 - fit_poly_log(f_leda, d2, 5), label=ant_ids[2])
    
    plt.xlim(40, 85)
    plt.minorticks_on()
    #plt.ylim(500, 7000)
    
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
