#!/usr/bin/env python
"""
# 09_fourier_fitting.py

Fourier/polynomial fitting of calibration parameters

"""
import matplotlib as mpl
import seaborn as sns
import tables as tb

from leda_cal.leda_cal import *
from leda_cal.skymodel import *
from leda_cal.dpflgr import *
from leda_cal.useful import poly_fit, fourier_fit


sns.set_style('ticks')
sns.set_context("paper",font_scale=1.5)

def quicklook(filename):
    
    balun_loss = hkl.load('cal_data/balun_loss.hkl')
    vna_cal    = hkl.load('cal_data/vna_calibration.hkl')
    
    ant_ids = ['a252x']#, 'a254x', 'a255x']
    
    plt.figure(figsize=(10, 4))
    for ant_id in ant_ids:
        ra = vna_cal[ant_id]["ra"]
        rl = vna_cal[ant_id]["rl"]
        f  = vna_cal["f"]
        F = compute_F(ra, rl)
        G = compute_G(rl)
        #G_f = fourier_fit(G, 0, 21)
        G_f = poly_fit(f, G, 11)
        F_f = fourier_fit(mag2(F), 0, 21)
        
        #ra_f = fourier_fit((1-mag2(ra)), 0, 21)
        ra_f = poly_fit(f, (1-mag2(ra)), 5)
        ra_f2 = fourier_fit((1-mag2(ra))-ra_f, 0, 21)
        
        ra_f = ra_f + ra_f2
        
        #plt.subplot(2, 3, 1)
        plt.subplot2grid((3, 3), (0, 0), rowspan=2)
        plt.plot(f, (1-mag2(ra)), c='#cc0000', label='$H_{\\rm{ant}}$')
        plt.plot(f, ra_f, c='#333333', label='$H_{\\rm{ant}}$')
        plt.yticks([0.5, 0.6, 0.7, 0.8])
        plt.ylim(0.5, 0.8)
        plt.xlim(40, 85)
        plt.xticks([40, 45, 50, 55, 60, 65, 70, 75, 80, 85], 
                   ["", "", "", "", "", "", "", "", "", ""])
        
        #plt.subplot(2,3,4)
        plt.subplot2grid((3, 3), (2, 0), rowspan=1)
        plt.plot(f, (1-mag2(ra))-ra_f, c='#333333')
        plt.xlim(40, 85)
        plt.xticks([40, 45, 50, 55, 60, 65, 70, 75, 80, 85],
                   [40, "", 50, "", 60, "", 70, "", 80, ""])
        plt.ylim(-0.002, 0.002)
        plt.yticks([-0.002, -0.001, 0, 0.001, 0.002])
        plt.xlabel("Frequency [MHz]")
        
        #plt.subplot(2, 3, 2)
        plt.subplot2grid((3, 3), (0, 1), rowspan=2)
        plt.plot(f, G, c='#cc0000', label='$H_{\\rm{lna}}$')
        plt.plot(f, G_f, c='#333333', label='$H_{\\rm{lna}}$')
        plt.yticks([0.997, 0.998, 0.999, 1.000])
        plt.ylim(0.997, 1.000)
        plt.xlim(40, 85)
        plt.xticks([40, 45, 50, 55, 60, 65, 70, 75, 80, 85],
                   ["", "", "", "", "", "", "", "", "", ""])
                   
        plt.subplot2grid((3, 3), (2, 1), rowspan=1)
        plt.plot(f, G-G_f, c='#333333')
        plt.xlim(40, 85)
        plt.xticks([40, 45, 50, 55, 60, 65, 70, 75, 80, 85],
                   [40, "", 50, "", 60, "", 70, "", 80, ""])
        plt.ylim(-0.0001, 0.0001)
        plt.xlabel("Frequency [MHz]")
                  
        #plt.subplot(2, 3, 3)
        plt.subplot2grid((3, 3), (0, 2), rowspan=2)
        plt.plot(f, mag2(F), c='#cc0000', label='$|F|^2$')
        plt.plot(f, F_f, c='#333333', label='$|F|^2$')
        plt.xlim(40, 85)
        plt.xticks([40, 45, 50, 55, 60, 65, 70, 75, 80, 85],
                   ["", "", "", "", "", "", "", "", "", ""])
        
                   
        #plt.subplot(2,3,6)
        plt.subplot2grid((3, 3), (2, 2), rowspan=1)
        plt.plot(f, mag2(F)-F_f, c='#333333')
        plt.xlim(40, 85)
        plt.xticks([40, 45, 50, 55, 60, 65, 70, 75, 80, 85],
                   [40, "", 50, "", 60, "", 70, "", 80, ""])
        plt.ylim(-0.002, 0.002)
        plt.yticks([-0.002, -0.001, 0, 0.001, 0.002])
        plt.xlabel("Frequency [MHz]")
    
    plt.tight_layout()
    #plt.legend(["$G$"], frameon=False, loc=4)
    #plt.text(40, 0.999, "$G$")
    plt.savefig("figures/cal-params.pdf")
    plt.show()

if __name__ == "__main__":
    
    import sys
    try:
        filename = sys.argv[1]
    except:
        print "USAGE: ./quicklook.py filename_of_hdf5_observation"
        exit()
    
    quicklook(filename)
