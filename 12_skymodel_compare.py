#!/usr/bin/env python
"""
# 12_skymodel_compare.py

Compare calibrated spectra against skymodel

"""
import matplotlib as mpl
import seaborn as sns
import tables as tb
from leda_cal.leda_cal import *
from leda_cal.skymodel import *
from leda_cal.useful import *
from leda_cal.dpflgr import *
from leda_cal.git import get_repo_fingerprint

sns.set_style('ticks')
sns.set_context("paper",font_scale=1.5)

def quicklook(filename, ant='252A', lfsm=False, emp=False, n_poly=7):
    h5 = tb.open_file(filename)

    T_ant = apply_calibration(h5)
    f_leda = T_ant['f']
    lst    = T_ant['lst']
    
    print T_ant.keys()
    
    ant_ids = [ant,]
      
    print("Plotting %s..." % ant_ids[0])
    fig, ax = plt.subplots(figsize=(8, 6))

    mid = closest(lst, 11)
    print T_ant['lst'][mid]
    sl  = 20
    
    T_flagged = rfi_flag(T_ant[ant_ids[0]], freqs=f_leda)
    
    sl  = 250
    
    plt.subplot2grid((3, 1), (0, 0), rowspan=2) 
    plt.plot(f_leda, T_flagged[mid-sl:mid+sl].mean(axis=0), c='#009933')
    
    
    import hickle as hkl
    gsm = hkl.load("cal_data/gsm-spec-lst11.hkl")
    plt.plot(gsm["f"], gsm["T_ew"], c='#333333', ls='dashed')
    
    if lfsm and emp:
        smdl = SkyModelLFSMEmp
        smlbl = 'LFSM+Emp'
    elif lfsm and not emp:
        smdl = SkyModelLFSM
        smlbl = 'LFSM'
    elif not lfsm and emp:
        smdl = SkyModelGSMEmp
        smlbl = 'GSM+Emp'
    else:
        smdl = SkyModelGSM        
        smlbl = 'GSM'
    s = smdl(pol='y' if ant_ids[0][-1] == 'A' else 'x')
    asm = s.generate_tsky(lst[mid-sl:mid+sl], f_leda*1e6).mean(axis=0)
    plt.plot(f_leda, asm)
    
    
    d = T_flagged[mid-sl:mid+sl].mean(axis=0)
    f_t, d_t = trim(f_leda, d, 40, 80)
    if ant_ids[0][-1] == 'A':
        T_hsm = np.interp(f_t, gsm["f"], gsm["T_ew"])
    else:
        T_hsm = np.interp(f_t, gsm["f"], gsm["T_ns"])
    T_asm = np.interp(f_t, f_leda, asm)
    
    scale_offset = np.mean(T_hsm / d_t)
    scale_offset_asm = np.mean(T_asm / d_t)
    print scale_offset, scale_offset_asm
    
    #plt.plot(f_)
    plt.xlim(40, 80)
    plt.ylim(500, 8000)
    plt.minorticks_on()
    plt.ylabel("Temperature [K]")
    plt.legend([ant_ids[0], "GSM (.hkl)", smlbl])
    
    plt.subplot2grid((3, 1), (2, 0), rowspan=1)
    plt.plot(rebin(f_t, 8), rebin(d_t / T_hsm, 8), c='#333333', linestyle='--')
    plt.plot(rebin(f_t, 8), rebin(scale_offset * d_t / T_hsm, 8), c='#333333', linestyle='--')
    plt.plot(rebin(f_t, 8), rebin(d_t / T_asm, 8), c='#333333')
    plt.plot(rebin(f_t, 8), rebin(scale_offset_asm * d_t / T_asm, 8), c='#333333')
    plt.xlabel("Frequency [MHz]")
    #plt.yticks([0.85, 0.87, 0.89, 0.91, 0.93])
    #  plt.ylim(0.85, 0.93)
    plt.ylabel("data / model")
    plt.tight_layout()
    plt.minorticks_on()
    plt.text(0.005, 0.005, get_repo_fingerprint(), transform=fig.transFigure, size=8)
    plt.savefig("figures/skymodel-compare.pdf")
    plt.show()
    
    resid = d_t - T_hsm
    resid -= fit_poly(f_t, resid, n_poly)
    resid_asm = d_t - T_asm
    resid_asm -= fit_poly(f_t, resid_asm, n_poly)
    
    plt.plot(f_t, resid, linestyle='--')
    plt.plot(f_t, resid_asm)
    plt.text(0.005, 0.005, get_repo_fingerprint(), transform=fig.transFigure, size=8)
    plt.show()

if __name__ == "__main__":
    import optparse, sys

    usage = '%prog [opts] filename_of_hdf5_observation'
    o = optparse.OptionParser()
    o.set_usage(usage)
    o.set_description(__doc__)
    o.add_option('--ant', dest='ant', action='store', default='252A', 
                 help='Name of the antenna to plot. Default: 252A')
    o.add_option('--lfsm', dest='lfsm', action='store_true', default=False,
                 help='Use the LFSM instead of the GSM')
    o.add_option('--empirical', dest='emp', action='store_true', default=False,
                 help='Apply an empirical corretion to the dipole gain pattern model')
    o.add_option('--n_poly', dest='n_poly', action='store', default=7, type="int",
                 help='Order of the polynomial to fit to the residuals. Default: 7')
    
    opts, args = o.parse_args(sys.argv[1:])
    
    if len(args) != 1:
        o.print_help()
        exit(1)
    else:
        filename = args[0]
        
    quicklook(filename, ant=opts.ant, lfsm=opts.lfsm, emp=opts.emp, n_poly=opts.n_poly)
