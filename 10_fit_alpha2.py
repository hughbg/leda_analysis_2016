#!/usr/bin/env python
"""
# 10_fit_alpha2.py

Fit sky power law (frequency^alpha) to the data over LST bins,
showing how it changes as a function of time.

NOTE: Also see 07_fit_alpha.py

"""
import matplotlib as mpl
import seaborn as sns
import tables as tb
from leda_cal.leda_cal import *
from leda_cal.skymodel import *
from leda_cal.dpflgr import *
from leda_cal.useful import fit_poly, closest, trim, fourier_fit, poly_fit, rebin

from lmfit import minimize, Parameters, report_fit

sns.set_style('ticks')
sns.set_context("paper",font_scale=1.5)

def residual(params, model, x, data, errs=1.0):
    model_vals = model(params, x)
    return (data - model_vals) / errs

def model(params, x):
    a0    = params['a0'].value
    b     = params['b'].value
    #A1     = params['A1'].value
    #A2     = params['A2'].value
    #A3     = params['A3'].value     
    #P2     = params['P2'].value

    model  = a0*np.power(x / 70, b) + 2.73 #+ A1 * x + A2 * x**2 + A3 * x**3 
    #print mo
    
    return model    

def quicklook(filename):
    
    n = 1
    fs, fe = 40, 83
    h5 = tb.open_file(filename)

    T_ant = apply_calibration(h5)
    f_leda = T_ant['f']
    lst_leda = T_ant["lst"]

    
    print T_ant.keys()
    
    ant_ids = ['252A']#, '254A', '255A', '254B', '255B']
      
    print("Plotting...")
    T0_all = []
    alpha_all = []
    lst_trials = np.arange(1, 23, 0.5)
    
    for ll in lst_trials:
    
        mid = closest(lst_leda, ll)
        sl  = 100
    
        
        for ant_id in ant_ids:
            T_flagged = T_ant[ant_id][mid-sl:mid+sl]
            print T_ant[ant_id].shape
            T_flagged = rfi_flag(T_flagged, thr_f=0.2, thr_t=0.2, rho=1.5,
                     bp_window_f=16, bp_window_t=16, 
                     max_frac_f=0.5, max_frac_t=0.5)    
            #T_flagged = np.ma.array(T_flagged)
    
            d = T_flagged.mean(axis=0)
            d_std = T_flagged.std(axis=0)

            print d.shape, f_leda.shape
        
            print ant_id, mid, sl

            f_t, d_t = trim(f_leda, d, fs, fe)
            f_t, d_std_t = trim(f_leda, d_std, fs, fe)
            f_t = np.ma.array(f_t)
            f_t.mask = d_t.mask
    
            params = Parameters()
            params.add('a0', value=1e6)
            params.add('b', value=-2.49)

            fc = f_t.compressed()
            dc = d_t.compressed()
            dc_errs = d_std_t.compressed()
        
            out = minimize(residual, params, args=(model, fc, dc, dc_errs))
    
            print report_fit(out)
            alpha_all.append(out.params['b'].value)
            T0_all.append(out.params['a0'].value)
            
            print "B: ", out.params['b'].value
    
            # print out.values
            #p = fit_poly(f_leda, d, n=n, print_fit=True)
    
            p = model(out.params, fc)
    
            print p.shape
            plt.figure("ANTZ")
            plt.plot(fc, dc, label=ant_id)
            #plt.plot(fc, p, label='fit')
    
            plt.xlim(fs, fe)
    plt.minorticks_on()
    plt.ylim(500, 7000)
    plt.xlabel("Frequency [MHz]")
    plt.ylabel("Temperature [K]")
    plt.legend(frameon=False)
    plt.show()
    

    alpha = np.array(alpha_all)
    print "mean std max min"
    print np.mean(alpha), np.std(alpha), np.max(alpha), np.min(alpha)
    plt.plot(lst_trials, alpha, c='#333333')
    plt.xlabel("Sidereal time [hr]")
    
    plt.ylabel("Spectral index $\\alpha$")
    plt.minorticks_on()
    plt.ylim(-2.4, -2.27)
    plt.twinx()
    plt.plot(lst_trials, T0_all, c='#cc0000', ls='dashed')
    plt.ylabel("$T_{70}$ [K]", color='#cc0000')
    plt.tight_layout()
    plt.minorticks_on()
    plt.xlim(0, 24)
    plt.ylim(1500, 3300)
    plt.savefig("figures/spectral-index.pdf")
    plt.show()
    
if __name__ == "__main__":
    
    import sys
    try:
        filename = sys.argv[1]
    except:
        print "USAGE: ./quicklook.py filename_of_hdf5_observation"
        exit()
    
    quicklook(filename)