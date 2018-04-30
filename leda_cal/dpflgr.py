#!/usr/bin/env python
"""
sumthreshold.py
---------------

First attempt at porting AOFlagger techniques to python.

TODO: Automated thresholding, better bandpass removal, etc.

Required modules:

* h5py - read data
* bottleneck - speedy moving averages
* matplotlib - plotting
"""

import datetime, os, sys
import h5py
import numpy as np
import pylab as plt
import bottleneck as bn

def fit_poly(x, y, n=5, log=False):
    """ Fit a polynomial to x, y data 
    
    x (np.array): x-axis of data (e.g. frequency)
    y (np.array): y-axis of data (e.g temperature)
    n (int): number of terms in polynomial (defaults to 5)
    """
    
    x_g = x
    x = np.ma.array(x, mask=y.mask).compressed()
    y = y.compressed()
    if log:
        yl = np.log10(y)
    else:
        yl = y
    fit = np.polyfit(x, yl, n)
    p = np.poly1d(fit)
    
    if log:
        return 10**(p(x_g))
    else:
        return p(x_g)


def plot_waterfall(d, freqs, lsts, t_unit='hr', f_unit='MHz',
                   raster=True):
    """ Plot imshow with LST and frequency on axis
    :param d: data
    :param freqs: frequency array
    :param lsts: LST array
    """

    plt.imshow(d, aspect='auto', interpolation='none',
               rasterized=raster,
               extent=(freqs[0], freqs[-1], lsts[0], lsts[-1]))
    plt.xlabel("Frequency [%s]" % f_unit)
    plt.ylabel("Time [%s]" % t_unit)

def check_mask(f):
    """ Decorator to check for masked arrays.
    
    Check that the first argument to a function is a masked array.
    If not, convert it into one.
    """
    def wrapper(*args, **kwargs):
        data = args[0]
        try:
            mask = data.mask
        except AttributeError:
            data = np.ma.array(data)
            mask = data.mask
            args = list(args)
            args[0] = data
            args = tuple(args)
        return f(*args, **kwargs)
    return wrapper

@check_mask
def simple_flag(data, threshold):
    """ Simple threshold flagging """
    
    d = np.diff(data)
    d = np.concatenate((np.array([0]), d))
    
    data.mask = np.abs(d) > threshold
    
    return data
           

@check_mask
def sum_threshold(data, thr_f, thr_t=None, scales=None, rho=1.5,
                  plot_progress=False, verbose=False):
    """ Apply Sum-Threshold method 
    
    This function applies a set ofmoving averages to the data along both 
    time and frequency axes, then checks if the output are above a threshold.
    This is the basic technique used in AOFlagger's algorithm.
    
    data (np.ma.array): data to flag, 2D array (time, freq)
    thr_f (int): threshold over which to flag on frequency axis
    thr_t (int): threshold over which to flag on time axis
    scales (list): list of window sizes (ints) to do moving average over.
            Defaults to None, in which case it uses [1,2,4,8,16,32,64] 
    rho (float): Threshold setting base. From eqn 12 in Offringa et. al. 2010:
                               thr_1 
                 thr_i =  --------------
                          rho^(log_2(i))
                 A value of 1.5 is suggested as being "empirically good"
    """

    if scales is None:
        scales = [1, 2, 4, 8, 16, 32, 64]
    
    if thr_t is None:
        thr_t = thr_f
    
    mask = np.copy(data.mask)
    
    thr1_f = thr_f
    thr1_t = thr_t
    
    # do first stage of flagging:
    mask_f = np.greater_equal(np.abs(data-1), thr_f)
    mask_t = np.greater_equal(np.abs(data-1), thr_t)
    #mask_b = np.greater_equal(np.abs(summed_b-1), np.sqrt(thr_f * thr_t))
    mask_s = np.logical_or(mask_f, mask_t)
    #mask_s = np.logical_or(mask_s, mask_b)        
    mask   = np.logical_or(data.mask, mask_s)
    data[mask] = np.sqrt(thr_f * thr_t)
    
    for window in scales:
        
        thr_f = thr1_f / np.power(rho, np.log2(window))
        thr_t = thr1_t / np.power(rho, np.log2(window))
        
        if window > 1:
            summed_f = bn.move_nanmean(data, window, axis=1)
            summed_t = bn.move_nanmean(data, window, axis=0)
            #summed_b = bn.move_nanmean(summed_f, int(np.sqrt(window)), axis=0)
        
            mask_f = np.greater_equal(np.abs(summed_f-1), thr_f)
            mask_t = np.greater_equal(np.abs(summed_t-1), thr_t)
            #mask_b = np.greater_equal(np.abs(summed_b-1), np.sqrt(thr_f * thr_t))
            mask_s = np.logical_or(mask_f, mask_t)
            #mask_s = np.logical_or(mask_s, mask_b)        
            mask   = np.logical_or(data.mask, mask_s)
            data[mask] = 1 + np.sqrt(thr_f * thr_t)
            data.mask  = mask
        else:
            summed_f = data
            summed_t = data

        if verbose:
            print "M: %i, Xi_f: %2.2e, Xi_t: %2.2e" % (window, thr_f, thr_t)

        if plot_progress:
            plt.figure()
            plt.subplot(221)
            plt.title("summed f: %i" % window)
            plt.imshow(summed_f, aspect='auto', interpolation='none', rasterized=True)
            plt.colorbar()
            plt.subplot(222)
            plt.title("summed t: %i" % window)
            plt.imshow(summed_t, aspect='auto', interpolation='none', rasterized=True)
            plt.colorbar()
            plt.subplot(223)
            plt.title("flagged: %i" % window)
            plt.imshow(data, aspect='auto', interpolation='none', rasterized=True)
            plt.colorbar()
    if plot_progress:
        plt.show()
        
    return data.mask

@check_mask
def estimate_bandpass(data, window_f=32, window_t=32):
    """ Estimate bandpass by rolling median over time
    
    data (np.ma.array): data array with axes (freq, time)
    window (int): size of moving window over which to compute
                  bandpass estimated by median.
    
    TODO: Fit a polynomial instead?
    """
        
    est = bn.move_median(data, window=window_f, axis=0)
    est = bn.move_median(est, window=window_t, axis=1)
    
    return est

@check_mask
def flag_absolute(data, thr_max=10, thr_min=0):
    """ Flag anything above or below absolute thresholds.
    
    First pass for horrific data.
    
    data (np.ma.array): data array with axes (freq, time)
    
    """
    data.mask = np.logical_or(data.mask, data > thr_max)
    data.mask = np.logical_or(data.mask, data < thr_min)
    return data.mask

@check_mask
def flag_fraction(data, max_frac_f=0.5, max_frac_t=0.5):
    """ Get rid of integs / freq channels with high occupancy
    
    data (np.ma.array): data array with axes (freq, time)
    max_frac_f (float): between 0 and 1, max occupancy fraction
                        for a given frequency channel
    max_frac_t (float): between 0 and 1, max occupancy fraction
                        for a given integration
    
    """
    occ_f = np.sum(data.mask, axis=0) / float(data.shape[0])
    occ_t = np.sum(data.mask, axis=1) / float(data.shape[1])
    
    bad_f = occ_f > max_frac_f
    bad_t = occ_t > max_frac_t
    
    data.mask[bad_t, :] = True
    data.mask[:, bad_f] = True
    
    return data.mask

@check_mask
def flag_window(data, window_f, window_t):
    """ Flag either side where stats aren't right
    
    The rolling window doesn't work properly on the edges, so flag this.
    
    data (np.ma.array): data array with axes (freq, time)
    window (int): size of window applied to data    
    
    """
    data.mask[:window_t, :] = True
    data.mask[-1 * window_t:, : ] = True
    data.mask[:, :window_f] = True
    data.mask[:, -1 * window_f:] = True
    return data.mask

def clip1(data, bp_window_f=8, bp_window_t=8):
  limit = 4
  for i in range(data.shape[1]):
    flat = bn.move_nanmean(data[:, i], bp_window_t)
    flat = data[:, i]-flat			# this will also insert the mask
    std = np.std(flat[np.logical_not(np.logical_or(np.isnan(flat), flat.mask))]); 
    clip_mask = np.logical_or(flat < -limit*std, limit*std < flat) 
    data.mask[:,i] = np.logical_or(data.mask[:,i], clip_mask)

  
def clip(data, bp_window_f=8, bp_window_t=8):

  # Get the standard deviation of the high (by frequency) third of the data, for clipping
  cut = data.shape[1]/3
  chunk = data[:, data.shape[1]-cut:]
  chunk = bn.move_nanmean(chunk, bp_window_t, axis=0)
  chunk = bn.move_nanmean(chunk, bp_window_f, axis=1)
  chunk = data[:, data.shape[1]-cut:]-chunk
  chunk = chunk[bp_window_t:, bp_window_f:]		# Because these edge values are nan now
  chunk = np.ravel(chunk)
  if np.ma.is_masked(chunk): chunk = chunk[chunk.mask==False]

  # Clipping values
  dmin = -4*np.std(chunk)
  dmax = 4*np.std(chunk); 

  # Mask the data. Have to flatten the data to find where to mask it
  flat = bn.move_nanmean(data, bp_window_t, axis=0)
  flat = bn.move_nanmean(flat, bp_window_f, axis=1)
  flat = data-flat;
  m = np.ma.mean(flat[bp_window_t:, bp_window_f:])			
  flat[:bp_window_t, :] = m		# Because these edge values are now Nan due to move_nanmean
  flat[:, :bp_window_f] = m
  flat -= m

  data.mask = np.logical_or(data.mask, flat>dmax)
  data.mask = np.logical_or(data.mask, flat<dmin)


@check_mask
def rfi_flag(data, thr_f, thr_t=None, rho=1.5,
             bp_window_f=64, bp_window_t=64, 
             max_frac_f=0.5, max_frac_t=0.5, scales=None, freqs=None):
    """ RFI Flagging routine for LEDA data
    
    Centered around a sum-threshold method, this method:
    1) Computes a moving estimate of the bandpass, and divides data through
       by this.
    2) Flags edges where moving estimate is dodgy
    3) Apply sum-threshold routine
    4) Flag any channel or integration with occupancy above a certain level.
    
    data (np.ma.array): data to flag, 2D array (time, freq)
    thr_max (float):   first stage flagging, maximum allowed value 
    thr_min (float):   first stage flagging, minimum allowed value 
    thr_f (float):     threshold over which to flag on frequency axis
    thr_t (float):     threshold over which to flag on time axis
    scales (list):   list of window sizes (ints) to do moving average over.
                     Defaults to None, in which case it uses [1,2,4,8,16,32,64] 
    med_window (int): size of window used to estimate bandpass
    max_frac_f (float): between 0 and 1, max occupancy fraction
                        for a given frequency channel
    max_frac_t (float): between 0 and 1, max occupancy fraction
                        for a given integration    
    """
    
    dtv_frequencies = [ 54.31, 60.31, 66.31, 76.31, 82.31 ]
    
    bpass = estimate_bandpass(data, window_f=bp_window_f, window_t=bp_window_t)
    to_flag = data / bpass
    #to_flag.mask = flag_absolute(to_flag, thr_max=thr_max, thr_min=thr_min)
    
    #plt.figure()
    #plt.imshow(data, aspect='auto', clim=(1,20))
    #plt.colorbar()
    #
    #plt.figure()
    #plt.imshow(bpass, aspect='auto')
    #plt.colorbar()
    #       
    #plt.figure()
    #plt.imshow(to_flag, aspect='auto')
    #plt.colorbar()
    #plt.show()
        
    to_flag.mask = sum_threshold(to_flag, thr_f, thr_t, scales)
    to_flag.mask = flag_fraction(to_flag, max_frac_f=max_frac_f, max_frac_t=max_frac_t)
    to_flag.mask = flag_window(to_flag, window_f=bp_window_f, window_t=bp_window_t)
    
    data.mask = to_flag.mask

  
    clip1(data, bp_window_f,bp_window_t)	# sigma clipping

    # DTV 
    if freqs is not None:
      for freq in dtv_frequencies:
        channel = len(freqs[freqs<freq])-1	# channel of the nearest freq
        channels = [ channel-1, channel, channel+1 ]
 
        for i in range(data.shape[0]):
          if np.ma.is_masked(data[i, channels]): 
            data.mask[i, channel-13:channel+250-13 ] = True	# The start of the TV subband is 13 channels back and 250 wide
        


    return data

    
