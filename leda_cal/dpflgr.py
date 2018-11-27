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

from filter import filter
from params import params

import robust

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
            data = np.ma.array(data, mask=np.zeros(data.shape, dtype=np.bool))
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
def sum_threshold(data, plot_progress=False, verbose=False):
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

    thr_f = params.thr_f
    thr_t = params.thr_t
    scales = params.scales
    rho = params.rho 
    
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
            summed_f = filter(data, window, axis=1)
            summed_t = filter(data, window, axis=0)
            #summed_b = filter(summed_f, int(np.sqrt(window)), axis=0, use_bn=use_bn)
        
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
def estimate_bandpass(data):
    """ Estimate bandpass by rolling median over time
    
    data (np.ma.array): data array with axes (freq, time)
    window (int): size of moving window over which to compute
                  bandpass estimated by median.
    
    TODO: Fit a polynomial instead?
    """
    
    est = filter(data, params.st_bp_window_f, axis=0)
    est = filter(est, params.st_bp_window_t, axis=1)
    
    return est

@check_mask
def flag_absolute(data):
    """ Flag anything above or below absolute thresholds.
    
    First pass for horrific data.
    
    data (np.ma.array): data array with axes (freq, time)
    
    """
    data.mask = np.logical_or(data.mask, data > params.thr_max)
    data.mask = np.logical_or(data.mask, data < params.thr_min)
    return data.mask

@check_mask
def flag_fraction(data):
    """ Get rid of integs / freq channels with high occupancy
    
    data (np.ma.array): data array with axes (freq, time)
    max_frac_f (float): between 0 and 1, max occupancy fraction
                        for a given frequency channel
    max_frac_t (float): between 0 and 1, max occupancy fraction
                        for a given integration
    
    """
    occ_f = np.sum(data.mask, axis=0) / float(data.shape[0])
    occ_t = np.sum(data.mask, axis=1) / float(data.shape[1])
    
    bad_f = occ_f > params.max_frac_f
    bad_t = occ_t > params.max_frac_t
    
    data.mask[bad_t, :] = True
    data.mask[:, bad_f] = True
    
    return data.mask

@check_mask
def flag_window(data):
    """ Flag either side where stats aren't right
    
    The rolling window doesn't work properly on the edges, so flag this.
    
    data (np.ma.array): data array with axes (freq, time)
    window (int): size of window applied to data    
    
    """
    data.mask[:params.st_bp_window_t, :] = True
    data.mask[-1 * params.st_bp_window_t:, : ] = True
    data.mask[:, :params.st_bp_window_f] = True
    data.mask[:, -1 * params.st_bp_window_f:] = True
    return data.mask

def clip(data):
    bp_window_t = params.sc_bp_window_t
    bp_window_f = params.sc_bp_window_f

    # Get the standard deviation of the high (by frequency) third of the data, for clipping
    cut = data.shape[1]/3
    chunk = data[:, data.shape[1]-cut:]
    chunk = bn.move_nanmean(chunk, bp_window_t, axis=0)
    chunk = np.roll(chunk, -bp_window_t/2+1, axis=0)
    chunk = bn.move_nanmean(chunk, bp_window_f, axis=1)
    chunk = np.roll(chunk, -bp_window_f/2+1, axis=1)
    chunk = data[:, data.shape[1]-cut:]-chunk
    chunk = chunk[bp_window_t:, bp_window_f:]		# Because these edge values are nan now
    chunk = np.ravel(chunk)
    if np.ma.is_masked(chunk):
        chunk = chunk[chunk.mask==False]
        
    # Clipping values
    dmin = -params.sigma*np.std(chunk)
    dmax = params.sigma*np.std(chunk)
    
    # Mask the data. Have to flatten the data to find where to mask it
    flat = bn.move_nanmean(data, bp_window_t, axis=0)
    flat = np.roll(flat, -bp_window_t/2+1, axis=0)
    flat = bn.move_nanmean(flat, bp_window_f, axis=1)
    flat = np.roll(flat, -bp_window_f/2+1, axis=1)
    flat = data-flat;
    m = np.ma.mean(flat[bp_window_t:, bp_window_f:])			
    flat[:bp_window_t, :] = m		# Because these edge values are now Nan due to move_nanmean
    flat[:, :bp_window_f] = m
    flat -= m
    
    data.mask = np.logical_or(data.mask, flat>dmax)
    data.mask = np.logical_or(data.mask, flat<dmin)


def clip1(data):
    nstart = np.ma.count(data)
    for i in xrange(data.shape[1]):
        flat = bn.move_nanmean(data[:, i], params.sc_bp_window_t, axis=-1)
        flat = np.roll(flat, -params.sc_bp_window_t/2+1, axis=-1)
        flat = data[:, i]-flat			# this will also insert the mask
        std = np.std(flat[np.logical_not(np.logical_or(np.isnan(flat), flat.mask))])
        if not np.isnan(float(std)):
            clip_mask = np.logical_or(flat < -params.sigma*std, params.sigma*std < flat) 
            data[:,i].mask = np.logical_or(data[:, i].mask, clip_mask)
            
    if np.ma.count(data) > nstart:
        raise RuntimeError("Number of points flagged went DOWN after clipping!")


@check_mask
def clip2(data, robust=True):
    """
    Alternate sigma-clipping flagger for spectrometer data.  This function
    assumes that the data have already been bandpassed in frequency and 
    time and then uses a iterative method to find and flag outliters.
    """
    
    for j in xrange(params.sc_passes):
        mask = data.mask*1
        
        for i in range(data.shape[1]):
            i0 = max([0, i-params.sc_bp_window_f/2])
            i1 = min([i+params.sc_bp_window_f/2, data.shape[1]-1])
            try:
                assert(robust)
                mn, st = robust.mean(data[:,i0:i1+1]), robust.std(data[:,i0:i1+1])
            except:
                mn, st = np.ma.mean(data[:,i0:i1+1]), np.ma.std(data[:,i0:i1+1])
            bad = np.where(np.abs(data[:,i]-1) > params.sigma*st)[0]
            mask[bad,i] |= True
                
        data.mask = mask*1
    return data.mask


def do_dtv_flagging(data, freqs):
    nstart = np.ma.count(data)
    for freq in params.dtv_frequencies:
        channel = len(freqs[freqs<freq])-1	# channel of the nearest freq
        channels = [ channel-1, channel, channel+1 ]
    
        for i in range(data.shape[0]):
            if np.ma.is_masked(data[i, channels]): 
                data[i, channel-13:channel+250-13 ].mask = True	# The start of the TV subband is 13 channels back and 250 wide
        
    if np.ma.count(data) > nstart:
        raise RuntimeError("Number of points flagged went DOWN after DTV flagging!")


@check_mask
def do_dtv_flagging2(data, freqs):
    """
    Alternate DTV flagger for spectrometer data.  This function uses the
    constrastin power between the DTV band edges (outer 0.25 MHz) and the 
    DTV band centers (inner 4.5 MHz) to identify DTV flares.
    """
    
    mask = data.mask*1
    dtv_times = []
 
    for ledge in (54, 60, 66, 76, 82):
         uedge = ledge + 6
         band = np.where( (freqs>=ledge) & (freqs<=uedge) )[0]
         trns = np.where( (freqs>=ledge+0.25) & (freqs<=uedge-0.25) )[0]
         empt = np.where( ((freqs>=ledge-0.25) & (freqs<ledge+0.25)) | ((freqs>uedge-0.25) & (freqs<=uedge+0.25)) )[0]
         
         pB = np.mean(data.data[:,band], axis=1)
         pT = np.mean(data.data[:,trns], axis=1)
         pE = np.mean(data.data[:,empt], axis=1)
         
         #import pylab
         #pylab.plot(pB-pE)
         #pylab.plot(pT-pE)
         #pylab.plot(pE-1)
         #pylab.plot(pE*0 + 3*pE.std())
         #pylab.show()
         
         st = np.std(pE)
         bad = np.where( np.abs(pT-pE) > 3*st )[0]
         for b in bad:
	     dtv_times.append(b)
             mask[b,band] |= True
             if b > 1:
                 mask[b-1,band] |= True
             if b < data.shape[0]-2:
                 mask[b+1,band] |= True
                 
    dtv_times = sorted(dtv_times)
    gap = -1
    where_gap = (-1, -1)
    for i in range(1, len(dtv_times)):
      this_gap = dtv_times[i]-dtv_times[i-1]
      if this_gap > gap:
        gap = this_gap
        where_gap = (dtv_times[i], dtv_times[i-1])

    data.mask = mask*1
    return data.mask, where_gap


@check_mask
def rfi_flag(data, freqs=None):
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
    
    
    if params.do_sum_threshold:
        try:
            to_flag
        except NameError:
            bpass = estimate_bandpass(data)
            to_flag = data / bpass
            
        to_flag.mask = sum_threshold(to_flag)
        to_flag.mask = flag_fraction(to_flag)
        to_flag.mask = flag_window(to_flag)
        
        data.mask = to_flag.mask
        
    if params.do_sigma_clip:
        try:
            to_flag
        except NameError:
            bpass = estimate_bandpass(data)
            to_flag = data / bpass
            
        to_flag.mask = clip2(to_flag)
        to_flag.mask = flag_fraction(to_flag)
        to_flag.mask = flag_window(to_flag)
        
        data.mask = to_flag.mask
        
    if params.do_dtv_flagging and freqs is not None:
        try:
            to_flag
        except NameError:
            bpass = estimate_bandpass(data)
            to_flag = data / bpass
            
        to_flag.mask, biggest_dtv_gap = do_dtv_flagging2(to_flag, freqs) 
        
        data.mask = to_flag.mask
        
    # Make sure all Nan are masked
    data.mask = np.logical_or(data.mask, np.isnan(data))
    
    return data, biggest_dtv_gap
