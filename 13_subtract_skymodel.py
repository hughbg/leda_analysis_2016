#!/usr/bin/env python
"""
# 13_subtract_skymodel.py

Subtract the skymodel from the calibrated data and plot.

NOTE: Old version, use 13b_subtract_skymodel.py

"""
import matplotlib as mpl
import seaborn as sns
import tables as tb
from leda_cal.leda_cal import *
from leda_cal.skymodel import *
from leda_cal.useful import *
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
    
    T_flagged = rfi_flag(T_ant[ant_ids[0]], thr_f=0.1, thr_t=0.1, rho=1.5,
             bp_window_f=8, bp_window_t=8, 
             max_frac_f=0.5, max_frac_t=0.5)
    
    sl  = 250
    
    n_poly = 7
    
    FA, FZ = 50, 85
    
    
    #plt.imshow(T_flagged, aspect='auto', cmap='viridis')
    #plt.show()
    
    #plt.subplot2grid((3, 1), (0, 0), rowspan=2) 
    #plt.plot(f_leda, T_flagged[mid-sl:mid+sl].mean(axis=0), c='#009933')
    
    
    import hickle as hkl
    gsm = hkl.load("gsm-spec-lst11.hkl")
    #plt.plot(gsm["f"], gsm["T_ew"], c='#333333', ls='dashed')
    
    d = T_flagged[mid-sl:mid+sl].mean(axis=0)
    
    d_flags = T_flagged.mask[mid-sl:mid+sl].sum(axis=0)
    
    d_errs  = T_flagged[mid-sl:mid+sl].std(axis=0) / np.sqrt(d_flags)
    
    #plt.figure("FLGS")
    #plt.plot(d_flags)
    #plt.show()
    
    f_t, d_t = trim(f_leda, d, FA, FZ)
    f_t, d_errs_t = trim(f_leda, d_errs, FA, FZ)
    T_ew = np.interp(f_t, gsm["f"], gsm["T_ew"])
    
    scale_offset = np.mean(T_ew / d_t)
    print scale_offset
    
    resid0 = scale_offset * d_t - T_ew
    
    print resid0
    
    model = poly_fit(f_t, resid0, n_poly, log=False)
    
    n_chan = 42
    A, B, C = rebin(f_t, n_chan), rebin(resid0, n_chan), rebin(d_errs_t, n_chan)
    print A.shape, B.shape, C.shape
    
    #print d.shape, d_flags.shape, d_errs.shape
    x = np.column_stack((A, B, C))
    np.savetxt("data_resids_gianni.txt", x)
    
    plt.figure("AAA")
    #plt.plot(f_t, resid0)
    #plt.plot(f_t, model)
    plt.plot(rebin(f_t, n_chan), rebin(resid0-model, n_chan))
    
    plt.figure("BBB")
    plt.plot(scale_offset * d_t - T_ew)
    plt.show()
    
    
    #plt.plot(f_)
    plt.xlim(FA, FZ)
    plt.ylim(500, 2000)
    plt.minorticks_on()
    plt.ylabel("Temperature [K]")
    plt.legend(["252A", "GSM"])
    
    plt.subplot2grid((3, 1), (2, 0), rowspan=1)
    plt.plot(rebin(f_t, n_chan), rebin(d_t / T_ew, n_chan), c='#333333')
    plt.plot(rebin(f_t, n_chan), rebin(scale_offset * d_t / T_ew, n_chan), c='#333333')
    plt.xlabel("Frequency [MHz]")
    #plt.yticks([0.85, 0.87, 0.89, 0.91, 0.93])
    #  plt.ylim(0.85, 0.93)
    plt.ylabel("data / model")
    plt.tight_layout()
    plt.minorticks_on()
    plt.savefig("figures/skymodel-compare.pdf")
    plt.show()
    
    resid = d_t - T_ew
    resid -= fit_poly(f_t, resid, n_poly)
    
    plt.plot(f_t, resid)
    plt.show()

if __name__ == "__main__":
    
    import sys
    try:
        filename = sys.argv[1]
    except:
        print "USAGE: ./quicklook.py filename_of_hdf5_observation"
        exit()
    
    quicklook(filename)