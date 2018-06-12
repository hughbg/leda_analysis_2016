"""
skymodel.py
-----------

Global Sky Model (GSM) class used for calibration of ledaspec data.

"""

import os
import hickle as hkl
import numpy as np
import pylab as plt

from scipy import interpolate
from scipy.optimize import curve_fit

class SkyModelGSM(object):
    """GSM sky model class using the NEC4 LWA dipole model"""
    
    # Relative path to the data
    _path = 'lfsm'
    
    # Common names for the .hkl files
    _name = 'gsm_nec'
    
    def __init__(self, pol='x', npol=7):
        self.pol = pol
        if pol == 'x':
            self.data = hkl.load(os.path.join(self._path, '%s_ew.hkl' % self._name))[:, 1:].T
            self.lsts = hkl.load(os.path.join(self._path, '%s_ew.hkl' % self._name))[:, 0].T

        if pol == 'y':
            self.data = hkl.load(os.path.join(self._path, '%s_ns.hkl' % self._name))[:, 1:].T
            self.lsts = hkl.load(os.path.join(self._path, '%s_ns.hkl' % self._name))[:, 0].T

        self.freqs = np.arange(30, 90.01, 5) * 1e6
        self.npol  = npol
        
        self.compute_gsm_splines()
    
    def curve_fit(self, f, d):
        """ Wrapper for curve_fit """
        f_l = np.log10(f / 70e6)
        d_l = np.log10(d)
        fit = np.polyfit(f_l, d_l, self.npol)
        
        return fit
    
    def compute_gsm_splines(self):
        """ Load GSM data and form interpolation splines """

        print "computing %s for Pol %s" % (self._name.upper(), self.pol)
        drift_data, drift_lsts, drift_freqs = self.data, self.lsts, self.freqs
        
        # Extend to full 24 hours then form interpolation spline
        nd = np.zeros((13, 145))
        nd[:, :144] = drift_data
        nd[:, 144]  = drift_data[:, 0]
        drift_lsts  = np.append(drift_lsts, drift_lsts[0]+24)
        drift_lsts[0] = 0.0
        drift_data  = nd
        
        fits = [ [] for ii in range(self.npol + 1)]
        #print fits
        for ii in range(len(drift_lsts)):

            fit = self.curve_fit(self.freqs, drift_data[:, ii])
            #if not ii%10:
                #print fit
            
            for jj in range(len(fit)):
                fits[jj].append(fit[jj])
        
        self.gsm_pols = np.array(fits)
        
        self.gsm_spline = []
        for kk in range(self.gsm_pols.shape[0]):
            self.gsm_spline.append(interpolate.interp1d(drift_lsts, self.gsm_pols[kk, :], kind='cubic'))

    def _generate_tsky(self, lst, freqs):
        try:
            poly = []
            for spl in self.gsm_spline:
                poly.append(spl(lst))

            p = np.poly1d(poly)
            return 10**(p(np.log10(freqs / 70e6)))
        except:
            return np.zeros_like(freqs)

    def generate_tsky(self, lst, freqs):
        """ Compute GSM Tsky for given lst and frequency range """

        if type(lst) is type(np.array([1,2])):
            tsky = np.zeros((len(lst), len(freqs)))
            ii = 0
            for l in lst:
                tsky[ii] = self._generate_tsky(l, freqs)
                ii += 1
            return tsky
        else:
            return self._generate_tsky(lst, freqs)

    def generate_tsky_vec(self, lst, freqs):
        """ TODO: Make this work. """
        tsky = np.zeros((len(lst), len(freqs)))
        np.apply_along_axis(self._generate_tsky, 0, tsky, lst, freqs)
        return tsky

class SkyModelGSMEmp(SkyModelGSM):
    """GSM sky model class the NEC4 LWA dipole model, plus empirical corrections"""
    
    # Relative path to the data
    _path = 'lfsm'
    
    # Common names for the .hkl files
    _name = 'gsm_emp'

class SkyModelLFSM(SkyModelGSM):
    """LFSM sky model class the NEC4 LWA dipole model"""
    
    # Relative path to the data
    _path = 'lfsm'
    
    # Common names for the .hkl files
    _name = 'lfsm_nec'

class SkyModelLFSMEmp(SkyModelGSM):
    """LFSM sky model class the NEC4 LWA dipole model, plus empirical corrections"""
    
    # Relative path to the data
    _path = 'lfsm'
    
    # Common names for the .hkl files
    _name = 'lfsm_emp'

# Backwards compatibility
SkyModel = SkyModelGSM

# Completeness
SkyModelGSMNEC = SkyModelGSM
SkyModelLFSMNEC = SkyModelLFSM

if __name__ == '__main__':        
    sgx = SkyModelGSM('x')
    sgy = SkyModelGSM('y')
    slx = SkyModelLFSM('x')
    sly = SkyModelLFSM('y')
    agx = SkyModelGSMEmp('x')
    agy = SkyModelGSMEmp('y')
    alx = SkyModelLFSMEmp('x')
    aly = SkyModelLFSMEmp('y')
    f = np.linspace(30, 90, 100) * 1e6

    gxx, gyy = [], []
    lxx, lyy = [], []
    axx, ayy = [], []
    bxx, byy = [], []

    for ii in range(0, 24):
        gxx.append(sgx.generate_tsky(ii, f))
        gyy.append(sgy.generate_tsky(ii, f))
        lxx.append(slx.generate_tsky(ii, f))
        lyy.append(sly.generate_tsky(ii, f))
        axx.append(agx.generate_tsky(ii, f))
        ayy.append(agy.generate_tsky(ii, f))
        bxx.append(alx.generate_tsky(ii, f))
        byy.append(aly.generate_tsky(ii, f))

    gxx, gyy = np.array(gxx), np.array(gyy)
    lxx, lyy = np.array(lxx), np.array(lyy)
    axx, ayy = np.array(axx), np.array(ayy)
    bxx, byy = np.array(bxx), np.array(byy)

    plt.plot(gxx[:, 0], linestyle='-', color='b')
    plt.plot(gyy[:, 0], linestyle='-', color='g')
    plt.plot(lxx[:, 0], linestyle='--', color='b')
    plt.plot(lyy[:, 0], linestyle='--', color='g')
    plt.plot(axx[:, 0], linestyle='-', color='r')
    plt.plot(ayy[:, 0], linestyle='-', color='k')
    plt.plot(bxx[:, 0], linestyle='--', color='r')
    plt.plot(byy[:, 0], linestyle='--', color='k')
    plt.show()