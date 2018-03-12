#!/usr/bin/env python

"""
plot_spa_all.py

Script for plotting SPA files from anritsu VNA
"""

import os
import numpy as np
from plot_spa import *
import skrf as rf
import glob

import h5py

def db20(x):
    return 20*np.log10(x)

def lin20(x):
    return 10**(x / 20.0)

def db10(x): 
    return 10*np.log10(x)
    
def lin10(x): 
    return 10.0**(x/ 10.0)

def y_factor(T_hot, T_cold, P_hot, P_cold):
    
    y = P_hot / P_cold
    T_dut = (T_hot - y * T_cold) / (y - 1)
    
    return T_dut

def enr_to_k(enr, T0=290, Tamb=290):
    """ Convert ENR to ENT 
    
    ENR = 10 Log10( (T - Tamb)/ T0 )
    
    """
    
    T = T0 * lin10(enr) + Tamb
    
    return T

def hp346a_enr(freqs, attenuation=0, cable_loss=-0.25, return_temps=False):
    """ Compute the ENR for a HP346A noise diode + attenuation
    
    This computes the ENR, taking into account cable losses and
    any extra attenuation added between the diode and the DUT.
    
    Parameters
    ----------
    freqs: np.array of frequency values, in MHz
    atten: np.array of attenuation values, in dB, 
           note: atten.shape == freq.shape
    cable_loss: np.array of cable losses, in dB
           note: cable_loss.shape == freq.shape
    return_temps: bool, return temperature (True) or ENR (False, default).
    
    Returns np.array of ENR in dB or equivalent temperature in K
    """
    f_cal = np.array([10, 100])
    enr_cal = np.array([5.37, 5.43])
    enr_fit = np.polyfit(f_cal, enr_cal, 1)
    enr     = np.poly1d(enr_fit)
    
    T_hp346 = enr_to_k(enr(freqs))
    
    L_atten = (1 - lin10(attenuation + cable_loss))
    T_atten = 290 * L_atten
    
    T_equiv = (1 - L_atten) * T_hp346 + T_atten
    
    #print L_atten.mean(), T_atten.mean(), T_equiv.mean()
    
    if not return_temps:
        enr_equiv   = db10(T_equiv / 290.0)
        return enr_equiv 
    else:
        return T_equiv

def hp346c_enr(freqs, attenuation=0, cable_loss=-0.25, return_temps=False):
    """ Compute the ENR for a HP346C noise diode + attenuation
    
    This computes the ENR, taking into account cable losses and
    any extra attenuation added between the diode and the DUT.
    
    Parameters
    ----------
    freqs: np.array of frequency values, in MHz
    atten: np.array of attenuation values, in dB, 
           note: atten.shape == freq.shape
    cable_loss: np.array of cable losses, in dB
           note: cable_loss.shape == freq.shape
    return_temps: bool, return temperature (True) or ENR (False, default).
    
    Returns np.array of ENR in dB or equivalent temperature in K
    """
    f_cal = np.array([10, 100])
    enr_cal = np.array([15.37, 15.43])
    enr_fit = np.polyfit(f_cal, enr_cal, 1)
    enr     = np.poly1d(enr_fit)
    
    T_hp346 = enr_to_k(enr(freqs))
    
    L_atten = (1 - lin10(attenuation + cable_loss))
    T_atten = 290 * L_atten
    
    #print T_atten
    
    T_equiv = (1 - L_atten) * T_hp346 + T_atten
    
    #print L_atten.mean(), T_atten.mean(), T_equiv.mean()
    
    if not return_temps:
        enr_equiv   = db(T_equiv / 290.0)
        return enr_equiv 
    else:
        return T_equiv


if __name__ == "__main__":
   
    
    f_hot       = 'y-factor/y-factor-y2-hot.spa'
    f_cold10    = 'y-factor/y-factor-y2-hot-10db.spa'
    f_cold6     = 'y-factor/y-factor-y2-hot-6db.spa'
    
    ff          = read_spa(f_hot)[0][:, 2]
    P_hot       = lin10(read_spa(f_hot)[0][:, 1])
    P_cold6     = lin10(read_spa(f_cold6)[0][:, 1])
    P_cold10    = lin10(read_spa(f_cold10)[0][:, 1])

    # Load attenuator and cable S21 measurements
    atten6 = rf.Network()
    atten6.read_touchstone('y-factor/yfactor-6db-pad.s2p')
    atten6_db = atten6.s21.s_db[:, 0, 0]
    atten6_f  = atten6.s21.f / 1e6    
    atten6_db = np.interp(ff, atten6_f, atten6_db)
    

    atten10 = rf.Network()
    atten10.read_touchstone('y-factor/yfactor-10db-pad.s2p')
    atten10_db = atten10.s21.s_db[:, 0, 0]
    atten10_f  = atten10.s21.f / 1e6
    atten10_db = np.interp(ff, atten10_f, atten10_db)  
    
    cable = rf.Network()
    cable.read_touchstone('y-factor/yfactor-cable-s21.s2p')  
    cable_db = cable.s21.s_db[:, 0, 0]
    cable_f  = cable.s21.f / 1e6
    cable_db = np.interp(ff, cable_f, cable_db)
    
    #plt.plot(ff, atten6_db)
    #plt.plot(ff, atten10_db)
    #plt.plot(ff, cable_db)
    #plt.show()
    
    T_hot    = hp346c_enr(ff, attenuation=0,          cable_loss=cable_db, return_temps=True)
    T_cold6  = hp346c_enr(ff, attenuation=atten6_db,  cable_loss=cable_db, return_temps=True)
    T_cold10 = hp346c_enr(ff, attenuation=atten10_db, cable_loss=cable_db, return_temps=True)
    
    plt.plot(ff, T_hot   , c='#339966', label='HP346C + 0dB pad')
    plt.plot(ff, T_cold6 , c='#663399', label='HP346C + 6dB pad')
    plt.plot(ff, T_cold10, c='#336699', label='HP346C + 10dB pad')
    plt.legend()
    plt.minorticks_on()
    plt.xlabel("Frequency [MHz]")
    plt.ylabel("Equivalent Noise Temperature [K]")
    plt.show()
    
    avg_temps = (np.average(T_hot), np.average(T_cold6), np.average(T_cold10))
    print "Hot: %2.2f K, 6dB: %2.2f K, 10db: %2.2f K" %  avg_temps
    
    plt.figure(figsize=(8, 6))                                                         
    plt.plot(ff, y_factor(T_hot, T_cold6, P_hot, P_cold6)      , c='#cc3333', label='Y1 0 - 6 dB')
    plt.plot(ff, y_factor(T_hot, T_cold10, P_hot, P_cold10)    , c='#666666', label='Y1 0 - 10 dB')
    plt.plot(ff, y_factor(T_cold6, T_cold10, P_cold6, P_cold10), c='#cc33cc', label='Y1 6 - 10 dB')

    f_hot       = 'y-factor/y-factor-y1-hot2.spa'
    f_cold10    = 'y-factor/y-factor-y1-hot2-10db.spa'
    f_cold6     = 'y-factor/y-factor-y1-hot2-6db.spa'

    P_hot       = lin10(read_spa(f_hot)[0][:, 1])
    P_cold6     = lin10(read_spa(f_cold6)[0][:, 1])
    P_cold10    = lin10(read_spa(f_cold10)[0][:, 1])
                                                           
    plt.plot(ff, y_factor(T_hot, T_cold6, P_hot, P_cold6)      , c='#339966', label='Y2 0 - 6 dB')
    plt.plot(ff, y_factor(T_hot, T_cold10, P_hot, P_cold10)    , c='#663399', label='Y2 0 - 10 dB')
    plt.plot(ff, y_factor(T_cold6, T_cold10, P_cold6, P_cold10), c='#336699', label='Y2 6 - 10 dB')

    plt.ylim(0, 5000)
    plt.xlim(25, 90)
    plt.minorticks_on()
    plt.legend()
    plt.xlabel("Frequency [MHz]")
    plt.ylabel("System temperature [K]")
    plt.savefig("boardy_rx_temp.png")


    plt.figure(figsize=(8, 6))
    plot_spa(f_hot,    show_fig=False, save_fig=False, label='hot',        c='#336699')
    plot_spa(f_cold6,  show_fig=False, save_fig=False, label='hot + 6dB',  c='#339966')
    plot_spa(f_cold10, show_fig=False, save_fig=False, label='hot + 10dB', c='#663399')
    
    plt.minorticks_on()
    plt.legend()
    
    plt.xlim(20, 95)
    plt.savefig("y_factor.png")
    plt.show()
