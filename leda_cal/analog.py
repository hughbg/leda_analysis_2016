#!/usr/bin/env python

import numpy as np
import pylab as plt
import hickle as hkl

def closest(xarr, val):
    return  np.argmin(np.abs(xarr - val))

def db(x): 
    return 10*np.log10(x)

def lin(x): 
    return 10**(x/20.0)

def polyfit(x, y, n=5):
    fit = np.poly1d(np.polyfit(x, y, n))
    return fit(x)

def s_to_t(ss):
    """ Convert S parameters into a T matrix"""
    pass 
    
def read_s2p(filename):
    """ Read an S2P file """
    
    ff = open(filename, 'r')
    dd = []
    
    # Read file
    lines = ff.readlines()
    
    # find last comment line (this is row headers)
    colnames_row = 0
    data_found = False
    for line in lines:
        if line.startswith('!') or line.startswith('#'):
            if not data_found:
                colnames_row += 1
            pass
        else:
            data_found = True
            dd.append(map(float, line.split()))
        
    colnames = np.array(lines[colnames_row-1].split())[1:]
    return colnames, np.array(dd)

def to_complex(data, linear=True, radians=False):
    """ Convert amp / phase data to complex."""
    r = data[:, 0]
    q = data[:, 1]
    
    if linear is True:
        r = 10**(r / 20)
    if radians is False:
        q = q / 180 * np.pi
    
    x = r * np.cos(q)
    y = r * np.sin(q)
    return x + 1j* y


def compute_F(ra, rl):
    """ Compute F-factor for Rogers calibration. 
    See Rogers and Bowman: Absolute Calibration.
    
    Parameters
    ----------
    ra: reflection coefficient of antenna
    rl: reflection coefficient of receiver
    """
    
    F = np.sqrt((1 - np.abs(rl)**2)) / (1 - ra*rl)
    return F

def compute_G(rl):
    """ Compute G-factor for Rogers calibration. 
    See Rogers and Bowman: Absolute Calibration.
    
    Parameters
    ----------
    rl: reflection coefficient of receiver
    """
    
    G = 1 - np.abs(rl)**2
    return G

def compute_T_ant(ra, rl, T_sky):
    """ Compute T_ant from eqn. 2 of Rogers and Bowman
    
    Parameters
    ----------
    ra: reflection coefficient of antenna
    rl: reflection coefficient of receiver   
    T_sky: sky temperature 
    """
    F = compute_F(ra, rl)
    T_ant = T_sky * (1 - np.abs(ra)**2) * np.abs(F)**2
    return T_ant
            
def compute_T_3p(ra, rl, T_sky):
    """ Compute T_3p from eqn. 13 of Rogers and Bowman
    
    Parameters
    ----------
    ra: reflection coefficient of antenna
    rl: reflection coefficient of receiver   
    T_sky: sky temperature 
    """
    F = compute_F(ra, rl)
    G = compute_G(rl)
    T_ant = T_sky * (1 - np.abs(ra)**2) * np.abs(F)**2 / G
    return T_ant

class AnalogComponent(object):
    """ S-parameters for a device under test. 
    
    Columns: Freq | S11 a/p | S21 a/p | S12 a/p | S22 a/p
    """
    def __init__(self, freq=None, s11=None, s21=None, s12=None, s22=None):
        self.freq = freq
        self.s11  = s11
        self.s21  = s21
        self.s12  = s12
        self.s22  = s22
    
    def read_s2p(self, filename):
        """ Read data in s2p file """
        colnames, dd = read_s2p(filename)
        
        self.freq = dd[:, 0]
        self.s11  = to_complex(dd[:, 1:3])
        self.s21  = to_complex(dd[:, 3:5])
        self.s12  = to_complex(dd[:, 5:7])
        self.s22  = to_complex(dd[:, 7:9])
    
    def _plot_mag_phase(self, s_param, c0='#cc0000', c1='#333333'):
        """ Plot s_param (helper fn)"""
        plt.subplot(211)
        plt.plot(self.freq, 20*np.log10(np.abs(s_param)), c=c0)
        plt.subplot(212)
        plt.plot(self.freq, np.rad2deg(np.angle(s_param)), c=c1)
    
    def plot_s11(self, show=True):
        plt.figure("S11")
        self._plot_mag_phase(self.s11)
        if show:
            plt.show()

    def plot_s12(self, show=True):
        plt.figure("S12")
        self._plot_mag_phase(self.s12)
        if show:
            plt.show()

    def plot_s21(self, show=True):
        plt.figure("S21")
        self._plot_mag_phase(self.s21)
        if show:
            plt.show()

    def plot_s22(self, show=True):
        plt.figure("S22")
        self._plot_mag_phase(self.s22)
        if show:
            plt.show()
    
    def plot_all(self, show=True):
        plt.figure("S-parameters", figsize=(12, 9))
        
        ii = 1
        for sid, sdata in [('S11', self.s11), ('S21', self.s21), ('S12', self.s12), ('S22', self.s22)]:
            plt.subplot(2, 2, ii)
            plt.title(sid)
            plt.plot(self.freq, 10*np.log10(np.abs(sdata)), c='#cc0000')
            plt.xlabel("Frequency [GHz]")
            plt.ylabel("Magnitude [dB]")
            plt.twinx()
            plt.plot(self.freq, np.rad2deg(np.angle(sdata)), c='#555555')
            plt.ylabel("Phase [deg]")
            plt.ylim(-185, 185)
            ii += 1
        
        plt.tight_layout()
        if show:
            plt.show()
        

def read_s2p_multifile(freq_fn, s11_fn, s21_fn, s12_fn, s22_fn):
    """ Combine forward and backward S-param measurements 
    
    Use to combine S11 and S21 and reverse S11 and S21, when
    VNA only has one port with reference input
    """
    af = AnalogComponent()
    af.read_s2p(freq_fn)

    a11 = AnalogComponent()
    a11.read_s2p(s11_fn)

    a21 = AnalogComponent()
    a21.read_s2p(s21_fn)

    a12 = AnalogComponent()
    a12.read_s2p(s12_fn)

    a22 = AnalogComponent()
    a22.read_s2p(s22_fn)
    
    a = AnalogComponent(af.freq, a11.s11, a21.s21, a12.s21, a22.s11)
    
    return a    