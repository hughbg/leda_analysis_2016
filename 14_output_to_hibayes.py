#!/usr/bin/env python
"""
# 14_output_to_hibayes.py

Outputs data to freqs.txt, spectrum.txt and spec_errs.txt for input into hibayes package
"""
import matplotlib as mpl
import seaborn as sns
import tables as tb
from leda_cal.leda_cal import *
from leda_cal.skymodel import *
from leda_cal.dpflgr import *
from leda_cal import useful

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
    
    n_chans = 42
    
    d0 = rfi_flag(d0, freqs=f_leda)
    
    d0 = useful.rebin(d0, 1, n_chans)
    f_leda = useful.rebin(f_leda, n_chans)
    
    d0_s = np.std(d0, axis=0) / np.sqrt(n_chans)**2
    d0 = np.median(d0, axis=0)
    
    

    #d1 = rfi_flag(d1, freqs=f_leda)
    #d1 = np.median(d1, axis=0)
    #
    #d2 = rfi_flag(d2, freqs=f_leda)
    #d2 = np.median(d2, axis=0)
    
    #plt.imshow(d0, cmap='viridis', aspect='auto')
    #plt.colorbar()
    #plt.show()
    #exit()
    
    return f_leda, d0, d0_s
    

if __name__ == "__main__":
    
    import sys
    try:
        filename = sys.argv[1]
    except:
        print "USAGE: ./quicklook.py filename_of_hdf5_observation"
        exit()
    
    f, d, d_s = quicklook(filename)
    
    np.savetxt("freqs.txt", f)
    np.savetxt("spectrum.txt", d)
    np.savetxt("spectrum_errs.txt", d_s)
    
    plt.subplot(2,1,1)
    plt.plot(f, d)
    plt.subplot(2,1,2)
    plt.plot(f, d_s)
    plt.show()
