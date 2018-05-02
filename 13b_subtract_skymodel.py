#!/usr/bin/env python
"""
# 13_subtract_skymodel.py

Subtract the skymodel from the calibrated data and plot.

NOTE: Never version, of 13_subtract_skymodel.py

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

F_START, F_STOP = 60, 82

def quicklook(filename):
    h5 = tb.open_file(filename)

    T_ant = apply_calibration(h5)
    f_leda = T_ant['f']
    lst    = T_ant['lst']
    
    print T_ant.keys()
    
    ant_ids, figout = ['252A', '254A', '255A'], "residuals-high-polA.pdf"
    #ant_ids, figout = ['BB', '254B', '255B'], "residuals-high-polB.pdf"
      
    print("Plotting...")
    fig, ax = plt.subplots(figsize=(8, 6))

    mid = closest(lst, 11)
    print T_ant['lst'][mid]
    
    
    import hickle as hkl
    gsm = hkl.load("cal_data/gsm-spec-lst11.hkl")
    sl  = 250    
    n_poly = 5
    n_chan = 42
    
    plt.figure("resid", figsize=(6, 8))

    for aa in ant_ids:
        if not aa == 'BB':
            T_flagged = rfi_flag(T_ant[aa], freqs=f_leda)

            #plt.plot(gsm["f"], gsm["T_ew"], c='#333333', ls='dashed')

            d = T_flagged[mid-sl:mid+sl].mean(axis=0)

            #d_flags = T_flagged.mask[mid-sl:mid+sl].sum(axis=0)
            #d_errs  = T_flagged[mid-sl:mid+sl].std(axis=0) / np.sqrt(d_flags)

            f_t, d_t = trim(f_leda, d, F_START, F_STOP)
            #f_t, d_errs_t = trim(f_leda, d_errs, 40, 80)
            T_ew = np.interp(f_t, gsm["f"], gsm["T_ew"])

            scale_offset = np.mean(T_ew / d_t)
            print scale_offset
            
            model_skypowerlaw = poly_fit(f_t, T_ew, 1, log=True)
            
            #resid0 = scale_offset * d_t - T_ew
            resid0 = d_t #scale_offset * d_t #- T_ew
            
            #plt.figure("m0")
            #plt.plot(T_ew - model_skypowerlaw)
            
            resid0 = resid0 #+ model_skypowerlaw
        
        pols = (1, 3, 5, 7)
        
        for ii, nn in enumerate(pols):
            plt.figure("resid")
            plt.subplot(len(pols), 1, ii+1)
            
            if aa == 'BB':
                plt.plot(0)
            else:
                if nn == 0:
                    model = 0
                else:
                    model = poly_fit(f_t, resid0, nn, log=True)
                    #model = fourier_fit(resid0, 0, nn)
                #model = 0

                #plt.plot(f_t, resid0)
                #plt.plot(f_t, model)
                
                plt.plot(rebin(f_t, n_chan), rebin(resid0-model, n_chan))
            plt.ylabel("Temperature [K]")
            plt.xlim(F_START, F_STOP)
            plt.minorticks_on()
            
    plt.xlabel("Frequency [MHz]")        
    
    plt.subplot(4,1,1,)
    plt.ylim(-100, 100)
    plt.subplot(4,1,2)
    plt.ylim(-30, 30)
    plt.subplot(4,1,3)
    plt.ylim(-30, 30)
    plt.subplot(4,1,4)
    plt.ylim(-12, 12)
    #plt.yticks([-40, -20, 0, 20, 40])
    
    plt.tight_layout()
    plt.savefig(figout)
    plt.show()

if __name__ == "__main__":
    
    import sys
    try:
        filename = sys.argv[1]
    except:
        print "USAGE: ./quicklook.py filename_of_hdf5_observation"
        exit()
    
    quicklook(filename)
