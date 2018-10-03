#!/usr/bin/env python
import scipy.interpolate
import tables as tb
import hickle as hkl
import ephem
from datetime import datetime
utcfromtimestamp = datetime.utcfromtimestamp
from analog import *
from useful import fourier_fit, poly_fit


def closest(xarr, val):
    return  np.argmin(np.abs(xarr - val))

def db(x): 
    return 10*np.log10(x)
    
def lin(x): 
    return 10**(x/10.0)

def polyfit(x, y, n=5):
    fit = np.poly1d(np.polyfit(x, y, n))
    return fit(x)

def mag2(x):
    """ Magnitude squared 
    
    Note magnitude is incorrectly stored in hickles with 10log10 def
    it should be 20log10 defn
    """
    return np.abs(x)**2

def extend(a, length):
    val = a[-1]
    num = length-len(a)
    x = np.append(a, np.zeros(num))
    x[-num:] = val
    return x

def interp_cplx(x_new, x_orig, d_orig):
    re = np.interp(x_new, x_orig, np.real(d_orig))
    im = np.interp(x_new, x_orig, np.imag(d_orig))
    d_new = re + 1j * im
    return d_new

def compute_F(ra, rl):
    F = np.sqrt((1 - mag2(rl))) / (1 - ra*rl)
    return F

def compute_H_lna(rl):
    G = 1 - mag2(rl)
    return G

def compute_GAV(ra, rl):
    G_AV = (1 - np.abs(ra)**2) / np.abs(1 - rl*ra)**2 * (1 - np.abs(rl)**2) 
    return G_AV

def compute_T_sky(ra, rl, T_3p, nw=0):
    """ Compute T_sky from T_3p calibrated spectra
    
    Args:
        ra (np.array): return loss of antenna
        rl (np.array): return loss of LNA
        nw (float): magnitude of noise wave
        T_3p (np.array): 3-state calibrated spectra (pre-VNA correction)
    """
    
    # Compute VNA cal factors
    H_lna   = compute_H_lna(rl)
    H_ant   = (1 - mag2(ra))
    F = compute_F(ra, rl)
    F_mag2  = mag2(F)
    
    # Fit models to them
    Fmag2_f  = fourier_fit(F_mag2, 0, 21)
    H_ant_f  = poly_fit(np.arange(H_ant.size), H_ant, 5)
    H_ant_f2 = fourier_fit(H_ant - H_ant_f, 0, 21)
    H_ant_f += H_ant_f2
    H_lna_f  = fourier_fit(H_lna, 0, 21)
    
    # Apply cal
    H_lna = extend(H_lna, T_3p.shape[1])
    H_ant = extend(H_ant, T_3p.shape[1])
    nw = extend(nw, T_3p.shape[1])
    F_mag2 = extend(F_mag2,  T_3p.shape[1])

    T_sky = (T_3p - nw) * H_lna / (H_ant * F_mag2) 
    return T_sky

def resample(f, data, new_f):
    if len(f) != len(data):
      print "Length of frequencies and data not the same in resample"
      exit(1)

    function = scipy.interpolate.interp1d(f, data)

    new_data = function(new_f)

    return new_data

def compute_noisewave(ant_id='a252y'):
    vna_cal    = hkl.load('cal_data/vna_calibration.hkl')
    
    ant_id = 'a252x'
    ra = vna_cal[ant_id]["ra"]
    rl = vna_cal[ant_id]["rl"]
    f  = vna_cal["f"]
    
    switch_s21 = hkl.load("cal_data/switch_s21.hkl")
    ra *= switch_s21
    
    F0 = compute_F(ra, rl)
    G = compute_G(rl)
    
    noisewave = hkl.load("cal_data/noisewave.hkl")
    T0 = noisewave[ant_id]['T0']
    Tu = noisewave[ant_id]['Tu']
    Tc = noisewave[ant_id]['Tc']
    Ts = noisewave[ant_id]['Ts']

    RA  = np.abs(ra)
    RA2 = RA**2 
    F2  = np.abs(F0)**2
    F   = np.abs(F0)
    PHI = np.angle(ra * F0)
    
    T_noise = T0 + Tu*RA2*F2 + RA*F*(Ts*np.sin(PHI) + Tc*np.cos(PHI))
    
    return T_noise

def apply_3ss_cal(ant_data, ant_id):
    T_rx_data = hkl.load('cal_data/rx_temperature_calibration.hkl')  # From ipython nb?
    try:
        T_h = T_rx_data[ant_id]['T_hot'][:-1]   # Check if -1 needed check all is right size
        T_c = T_rx_data[ant_id]['T_cold'][:-1]    
        T_c = extend(T_c, ant_data.shape[1])
        T_h = extend(T_h, ant_data.shape[1])
        T_ant = (T_h - T_c) * ant_data + T_c
    except ValueError:
        T_h = T_rx_data[ant_id]['T_hot']
        T_c = T_rx_data[ant_id]['T_cold']
        T_ant = (T_h - T_c) * ant_data + T_c        
    
    return T_ant

def apply_vna_cal(T_3p, ant_id):
    balun_loss = hkl.load('cal_data/balun_loss.hkl')    # Can plot this and see smooth line
    vna_cal    = hkl.load('cal_data/vna_calibration.hkl')
    
    ra = vna_cal[ant_id]["ra"][:-1]
    rl = vna_cal[ant_id]["rl"][:-1]
    L  = 10**(-balun_loss[ant_id] / 20.0)
    
    #T_3p = T_3p[:, :2290]
    #L = L[:2290]
    
    nw = compute_noisewave()
    T_sky_meas = compute_T_sky(ra, rl, T_3p, nw)
    
    L = extend(L, T_sky_meas.shape[1])
    T_sky_meas_corr = (T_sky_meas - 290 * (1 - L)) / L
  
    return T_sky_meas_corr


def apply_calibration(h5, apply_3ss=True, apply_vna=True):
    """ Apply 3SS calibration (convert to temperature) and apply VNA calibration """
    
    (latitude, longitude, elevation) = ('37.2397808', '-118.2816819', 1183.4839)

    ov = ephem.Observer()
    ov.lon = longitude
    ov.lat = latitude
    ov.elev = elevation
            
    tstamps = h5.root.data.cols.timestamp[:]
    f_leda = hkl.load('cal_data/vna_calibration.hkl')["f"]
    if f_leda.shape[0] > h5.root.data.cols.ant252_x.shape[1]:
        f_leda = f_leda[:-1]
    elif f_leda.shape[0] < h5.root.data.cols.ant252_x.shape[1]:
        last_index = len(f_leda)-1
        last_val = f_leda[-1]
        f_leda = extend(f_leda, h5.root.data.cols.ant252_x.shape[1])
        for i in range(last_index+1, len(f_leda)):
          f_leda[i] = last_val+(i-last_index)*.024
 
    
    lst_stamps = np.zeros_like(tstamps)
    utc_stamps = []
    for ii, tt in enumerate(tstamps):
        utc = utcfromtimestamp(tt)
        ov.date = utc
        lst_stamps[ii] = ov.sidereal_time() * 12.0 / np.pi
        utc_stamps.append(utc)
    
    if apply_3ss:
        print("Applying 3-state switching calibration")
        a252x = apply_3ss_cal(h5.root.data.cols.ant252_x[:], 'a252x')
        a254x = apply_3ss_cal(h5.root.data.cols.ant254_x[:], 'a254x')
        a255x = apply_3ss_cal(h5.root.data.cols.ant255_x[:], 'a255x')
        a252y = apply_3ss_cal(h5.root.data.cols.ant252_y[:], 'a252y')
        a254y = apply_3ss_cal(h5.root.data.cols.ant254_y[:], 'a254y')
        a255y = apply_3ss_cal(h5.root.data.cols.ant255_y[:], 'a255y')
    else:
        a252x = h5.root.data.cols.ant252_x[:]
        a254x = h5.root.data.cols.ant254_x[:]
        a255x = h5.root.data.cols.ant255_x[:]
        a252y = h5.root.data.cols.ant252_y[:]
        a254y = h5.root.data.cols.ant254_y[:]
        a255y = h5.root.data.cols.ant255_y[:]
    if apply_vna:
        print("Applying VNA calibration")
        a252x = apply_vna_cal(a252x, 'a252x')
        a254x = apply_vna_cal(a254x, 'a254x')
        a255x = apply_vna_cal(a255x, 'a255x')
        a252y = apply_vna_cal(a252y, 'a252y')
        a254y = apply_vna_cal(a254y, 'a254y')
        a255y = apply_vna_cal(a255y, 'a255y')
    
    print("Sorting by LST")
    sort_idx = lst_stamps.argsort()
    lst_stamps = lst_stamps[sort_idx]
    a252x = a252x[sort_idx]
    a254x = a254x[sort_idx]
    a255x = a255x[sort_idx]
    a252y = a252y[sort_idx]
    a254y = a254y[sort_idx]
    a255y = a255y[sort_idx]    
    utc_stamps = [ utc_stamps[sort_idx[i]] for i in range(len(sort_idx)) ]

       
    ants = {'lst': lst_stamps, 'utc' : utc_stamps, 'f' : f_leda,
            '252A': a252x, '252B': a252y,
            '254A': a254x, '254B': a254y,
            '255A': a255x, '255B': a255y,
            }
            
    return ants

def calibrate(data, caldata):
    T_H = caldata['T_H']
    T_C = caldata['T_C']
    G_S = caldata['G_S']
    T_NW = caldata['T_NW']
    # S   = caldata['scale']
    # O   = caldata['offset']

    # D = S * ((T_H - T_C) * data + T_C) / G_S + O
    D = ((T_H - T_C) * data + T_C) / G_S
    return D




def apply_new_calibration(h5):
    print "Using 2018 calibration method"

    (latitude, longitude, elevation) = ('37.2397808', '-118.2816819', 1183.4839)

    ov = ephem.Observer()
    ov.lon = longitude
    ov.lat = latitude
    ov.elev = elevation


    print "Loading cal files"
    h_252x = hkl.load('cal_data/cal_data_252A_fit.hkl')		# New cal files generated by leda_cal2 repo
    h_254x = hkl.load('cal_data/cal_data_254A_fit.hkl')
    h_255x = hkl.load('cal_data/cal_data_255A_fit.hkl')
    h_252y = np.zeros_like(h_255x)
    h_254y = hkl.load('cal_data/cal_data_254B_fit.hkl')
    h_255y = hkl.load('cal_data/cal_data_255B_fit.hkl')

    if not ( len(h_254x['f_mhz']) == len(h_252x['f_mhz']) and len(h_255x['f_mhz']) == len(h_252x['f_mhz']) and len(h_254y['f_mhz']) == len(h_252x['f_mhz']) 
		and len(h_255y['f_mhz']) == len(h_252x['f_mhz']) ):
	print "Consistency error: cal files have different frequencies"
	exit(1)

    f_leda = h_254x['f_mhz']

    tstamps = h5.root.data.cols.timestamp[:]

    if f_leda.shape[0] > h5.root.data.cols.ant252_x.shape[1]:
      f_leda = f_leda[:-1]
    elif f_leda.shape[0] < h5.root.data.cols.ant252_x.shape[1]:
      last_index = len(f_leda)-1
      last_val = f_leda[-1]
      print "Extending data to match wider frequency range"
      f_leda = extend(f_leda, h5.root.data.cols.ant252_x.shape[1])
      for i in range(last_index+1, len(f_leda)):
        f_leda[i] = last_val+(i-last_index)*.024
 
    lst_stamps = np.zeros_like(tstamps)
    utc_stamps = []
    for ii, tt in enumerate(tstamps):
        utc = utcfromtimestamp(tt)
        ov.date = utc
        lst_stamps[ii] = ov.sidereal_time() * 12.0 / np.pi
        utc_stamps.append(utc)

    print "Applying cal files"
    a252x = calibrate(h5.root.data.cols.ant252_x[:], h_252x)
    a252y = np.zeros_like(a252x)
    a254x = calibrate(h5.root.data.cols.ant254_x[:], h_254x)
    a255x = calibrate(h5.root.data.cols.ant255_x[:], h_255x)
    a254y = calibrate(h5.root.data.cols.ant254_y[:], h_254y)
    a255y = calibrate(h5.root.data.cols.ant255_y[:], h_255y)
    
    print("Sorting by LST")
    sort_idx = lst_stamps.argsort()
    lst_stamps = lst_stamps[sort_idx]
    a252x = a252x[sort_idx]
    a254x = a254x[sort_idx]
    a255x = a255x[sort_idx]
    a252y = a252y[sort_idx]
    a254y = a254y[sort_idx]
    a255y = a255y[sort_idx]    
    utc_stamps = [ utc_stamps[sort_idx[i]] for i in range(len(sort_idx)) ]

    ants = {'lst': lst_stamps, 'utc' : utc_stamps, 'f' : f_leda,
            '252A': a252x, '252B': a252y,
            '254A': a254x, '254B': a254y,
            '255A': a255x, '255B': a255y,
            }
            
    return ants

