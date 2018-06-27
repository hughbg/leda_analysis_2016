#!/usr/bin/env python
"""
# 07_fit_alpha.py

Fit sky power law (frequency^alpha) to the data over 10 MHz chunks,
starting at 40-50 MHz and ending at 70-80 MHz, then plot alpha as
a function of frequency.

NOTE: Also see 10_fit_alpha2.py

"""
import matplotlib as mpl
import seaborn as sns
import tables as tb
from leda_cal.leda_cal import *
from leda_cal.skymodel import *
from leda_cal.dpflgr import *
from leda_cal.useful import fit_poly, closest, trim, fourier_fit, poly_fit, rebin
from leda_cal.git import get_repo_fingerprint

from lmfit import minimize, Parameters, report_fit
sns.set_style('ticks')
sns.set_context("paper",font_scale=1.5)

def residual(params, model, x, data, errs=1.0):
    model_vals = model(params, x)
    return (data - model_vals) / errs

def model(params, x):
    amp    = params['amp'].value
    alpha  = params['alpha'].value
    
    model = amp*np.power(x, alpha) #+ off #+ a_s * np.sin(theta * x) #+ a_c * np.cos(theta * x)
    #print mo
    
    return model    

def quicklook(filename):
    
    n = 1
    fs, fe = 40, 80
    h5 = tb.open_file(filename)

    T_ant = apply_calibration(h5)
    f_leda = T_ant['f']
    
    print T_ant.keys()
    
    ant_ids = ['252A', '254A', '255A']
    pol_id  = 'y'
      
    print("Plotting...")
    fig, ax = plt.subplots(figsize=(8, 6))

    mid = T_ant["252A"].shape[0]/2
    
    print T_ant['lst'][mid]
    sl  = 250

    f0 = np.arange(40, 70)
    f1 = f0 + 10
    f_trials = zip(f0, f1)
    f_trials.append((40, 80))
    alpha_out = []

    for fs, fe in f_trials:
        T_flagged = T_ant[ant_ids[0]][mid-sl:mid+sl]
        T_flagged = rfi_flag(T_flagged, freqs=f_leda)

        d      = np.median(T_flagged, axis=0)
        d_errs = np.std(T_flagged, axis=0)

        f_trim, d = trim(f_leda, d, fs, fe)
        f_trim = np.ma.array(f_trim)
        f_trim.mask = d.mask

        f_trim2, d_errs = trim(f_leda, d_errs, fs, fe)
        d_errs = np.ma.array(d_errs)
        d_errs.mask = d.mask

        params = Parameters()
        params.add('amp', value=1e8, min=1e4, max=1e10)
        params.add('alpha', value=-2.5, min=-4, max=-1)

        fc = f_trim.compressed()
        dc = d.compressed()
        dc_errs = d_errs.compressed()

        print dc.shape, dc_errs.shape
        out = minimize(residual, params, args=(model, fc, dc, dc_errs))

        print report_fit(out)

        # print out.values
        #p = fit_poly(f_leda, d, n=n, print_fit=True)

        p = model(out.params, fc)
        alpha_out.append(out.params['alpha'].value)

        """
        plt.plot(fc, dc, label=ant_ids[0])
        plt.plot(fc, p, label='fit')

        plt.xlim(fs, fe)
        plt.minorticks_on()
        plt.ylim(500, 7000)

        ax.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
        ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())

        plt.xlabel("Frequency [MHz]")
        plt.ylabel("Temperature [K]")

        plt.legend(frameon=False)
        plt.show()
        """
    print "HERE", f0.shape, len(alpha_out)
    plt.plot(f0, alpha_out[:-1])

    plt.figure()
    res = dc - p
    res_1p = dc - poly_fit(fc, dc, 1, log=True)
    #res_2p = dc - poly_fit(fc, dc, 2, log=True)
    res_3p = dc - poly_fit(fc, dc, 3, log=True)
    #res_4p = dc - poly_fit(fc, dc, 4, log=True)
    res_5p = dc - poly_fit(fc, dc, 5, log=True)
    res_7p = dc - poly_fit(fc, dc, 7, log=True)
    #res_3p3f = res_3p - fourier_fit(res_3p, 0, 3)
    
    #plt.plot(fc, res)
    plt.plot(rebin(fc, 8), rebin(res_1p, 8))
    plt.plot(rebin(fc, 8), rebin(res_3p, 8))
    plt.plot(rebin(fc, 8), rebin(res_5p, 8))
    plt.plot(rebin(fc, 8), rebin(res_7p, 8))
    #plt.plot(rebin(fc, 8), rebin(res_5p, 8))
    
    print np.std(rebin(res_1p, 8))
    print np.std(rebin(res_3p, 8))
    print np.std(rebin(res_5p, 8))
    print np.std(rebin(res_7p, 8))
    #print np.std(rebin(res_5p, 8))
    #plt.plot(rebin(fc, 8), rebin(res_3p3f, 8))
    
    plt.text(0.005, 0.005, get_repo_fingerprint(), transform=fig.transFigure, size=8)
    
    plt.show()


if __name__ == "__main__":
    
    import sys
    try:
        filename = sys.argv[1]
    except:
        print "USAGE: ./quicklook.py filename_of_hdf5_observation"
        exit()
    
    quicklook(filename)
