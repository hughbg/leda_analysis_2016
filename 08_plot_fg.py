#!/usr/bin/env python
"""
# 08_plot_fg.py

Plot the calibration coefficients from Bowman et al (2012), F and G.

These coefficients are used when applying VNA calibration.

"""
import matplotlib as mpl
import seaborn as sns
import tables as tb
from leda_cal.leda_cal import *
from leda_cal.skymodel import *
from leda_cal.dpflgr import *

from lmfit import minimize, Parameters, report_fit

sns.set_style('ticks')
sns.set_context("paper",font_scale=1.5)



def quicklook():
    
    balun_loss = hkl.load('cal_data/balun_loss.hkl')
    vna_cal    = hkl.load('cal_data/vna_calibration.hkl')
    
    ant_ids = ['a252x']#, 'a254x', 'a255x']
    for ant_id in ant_ids:
        ra = vna_cal[ant_id]["ra"]
        rl = vna_cal[ant_id]["rl"]
        f  = vna_cal["f"]
        F = compute_F(ra, rl)
        G = compute_G(rl)
        plt.subplot(3,1,1)
        plt.plot(f, (1-mag2(ra)), c='#002147', label='$H_{\\rm{ant}}$')
        plt.yticks([0.5, 0.6, 0.7, 0.8])
        plt.ylim(0.5, 0.8)
        plt.xticks([40, 45, 50, 55, 60, 65, 70, 75, 80, 85], 
                   ["", "", "", "", "", "", "", "", "", ""])
        plt.subplot(3, 1, 3)
        plt.plot(f, G, c='#002147', label='$H_{\\rm{lna}}$')
        plt.yticks([0.997, 0.998, 0.999, 1.000])
        plt.ylim(0.997, 1.000)
        plt.subplot(3, 1, 2)
        plt.plot(f, mag2(F), c='#002147', label='$|F|^2$')
        plt.xticks([40, 45, 50, 55, 60, 65, 70, 75, 80, 85], 
                   ["", "", "", "", "", "", "", "", "", ""])

    
    for ii in (1,2,3):
        plt.subplot(3,1,ii)
        plt.xlim(40, 85)
        plt.rcParams["legend.fontsize"] = 12
        plt.legend(loc=4, frameon=False)
    #plt.legend(["$|F|^2$"], frameon=False, loc=4)
    
    plt.xlabel("Frequency [MHz]")
    plt.tight_layout()
    #plt.legend(["$G$"], frameon=False, loc=4)
    #plt.text(40, 0.999, "$G$")
    plt.show()

if __name__ == "__main__":

    quicklook()