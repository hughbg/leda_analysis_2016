#!/usr/bin/env python
"""
# plot_ds9.py

Plot data in a FITS format using DS9.

"""

import os
import seaborn as sns
from leda_cal.leda_cal import *
from leda_cal.dpflgr import *
import scipy.stats
import scipy.ndimage.filters
from leda_cal.useful import statistics, ensure_mask
from leda_cal.params import params
from leda_cal import lst_timing


fits_header = "SIMPLE  =                    T / file does conform to FITS standard             BITPIX  =                  -32 / number of bits per data pixel                  NAXIS   =                    2 / number of data axes                            NAXIS1  =                 XXXX / length of data axis 1                          NAXIS2  =                 YYYY / length of data axis 2                          EXTEND  =                    T / FITS dataset may contain extensions            COMMENT   FITS (Flexible Image Transport System) format is defined in 'AstronomyCOMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H CTYPE2  = 'LST     '                                                            CTYPE1  = 'FREQ    '                                                            CRPIX1  =                   1.                                                  CRPIX1  =                   1.                                                  CRVAL2  =               STARTL / LST                                            CRVAL1  =               STARTF / Frequency                                      CDELT2  =           LLLLLLLLLL /                                                CDELT1  =           FFFFFFFFFF /                                                END                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "

fits_header += " "*2880
fits_header = fits_header[:2880]

ovro_location = ('37.2397808', '-118.2816819', 1183.4839)
ovro = ephem.Observer(); (ovro.lat, ovro.lon, ovro.elev) = ovro_location
sun = ephem.Sun()

gal_center = ephem.FixedBody()  
gal_center._ra  =  '17 45 40.04'
gal_center._dec = '-29 00 28.1'
gal_center.name = "Galactic Center"

REQUIRED = 2400		# channels. Pad if necessary.

def rms_filter(data):
    return np.std(data[np.logical_not(np.isnan(data))])  

def biggest_gap(times):	# times must be sorted
    gap = -1
    where_gap = (-1, -1)
    for i in range(1, len(times)):
      this_gap = times[i]-times[i-1]
      if this_gap > gap:
        gap = this_gap
        where_gap = (times[i], times[i-1])
    return where_gap

def pad_data(data):

  if data.shape[1] == REQUIRED: return data

  missing = REQUIRED-data.shape[1]
  if missing > 0:
    pad = np.ma.masked_all((data.shape[0], missing))
    new_data = np.ma.concatenate((ensure_mask(data), pad), axis=1)

  return new_data



def quicklook(filename, ant, flag, noise, no_show, all_lsts, new_cal):
    global fits_header
    
    h5 = tb.open_file(filename)
    if new_cal: T_ant = apply_new_calibration(h5)
    else: T_ant = apply_calibration(h5)
    f_leda = T_ant['f']

    lst_stamps = T_ant['lst']
    utc_stamps = T_ant['utc']

    print "Start", lst_stamps[0], utc_stamps[0], lst_stamps[-1], utc_stamps[-1]

    
    # Report discontinuities in time
    for i in range(1,len(lst_stamps)):
        if lst_stamps[i]-lst_stamps[i-1] > 1/60.0:	# 1 minute
            print "Discontinuity at LST", lst_stamps[i], (lst_stamps[i]-lst_stamps[i-1])*60*60, "seconds"
            
   # Delete sun up LSTS
    if not all_lsts:
       time_info = lst_timing.LST_Timing(lst_stamps, utc_stamps)
       border_bottom, night_bottom, night_top, border_top = time_info.calc_night()
       if not border_top:
	 raise RuntimeError("No LSTs available at night time (use --all_lsts to see all)")
       print "Cutting out Sun/Galaxy up. LSTS in DS9 may be WRONG!!!!!!"
       lst_stamps = lst_stamps[night_bottom:night_top]
       utc_stamps = lst_stamps[night_bottom:night_top]
       print len(lst_stamps), "usable LSTs"
    else:
        print "Using all LSTs"

    if len(lst_stamps) == 0:
        raise RuntimeError("There is no data to display (number of LSTs is 0)")
        
    fits_header = fits_header.replace("XXXX", ( "%4d" % REQUIRED ))
    fits_header = fits_header.replace("YYYY", ( "%4d" % len(lst_stamps) ))

    lst_diff = (str((lst_stamps[0]-lst_stamps[-1])/(len(lst_stamps)-1))+"0000000000")[:10]
    freq_diff = (str((f_leda[-1]-f_leda[0])/(len(f_leda)-1))+"0000000000")[:10]
    fits_header = fits_header.replace("LLLLLLLLLL", lst_diff)
    fits_header = fits_header.replace("FFFFFFFFFF", freq_diff)

    lst_start = (str(lst_stamps[-1])+"000000")[:6]			# start and the end for FITs in our orientation
    freq_start = (str(f_leda[0])+"000000")[:6]
    fits_header = fits_header.replace("STARTL", lst_start)
    fits_header = fits_header.replace("STARTF", freq_start)
    
    if flag:
        print "Flagging"
        if not all_lsts:
	  data = rfi_flag(T_ant[ant][border_bottom:border_top], freqs=f_leda)[0]
	  data = data[night_bottom-border_bottom:night_top-border_bottom]	# Remove border
          total = data.shape[0]*data.shape[1]
          num_in = np.ma.MaskedArray.count(data)
          print ant, ( "%.1f%%" % (100*float(total-num_in)/total) ), "flagged.", "Count:", total-num_in

	else:
          data, dtv_tms = rfi_flag(T_ant[ant], freqs=f_leda)
	  print ant, "Biggest DTV gap", lst_stamps[biggest_gap(dtv_tms)[1]], "-", lst_stamps[biggest_gap(dtv_tms)[0]]

    else:
        if not all_lsts:
          data = T_ant[ant][night_bottom:night_top]
        else:
          data = T_ant[ant]

            
    #statistics(data)
    
    # Test the image is oriented right
    #for i in range(200):
    #  for j in range(200): data[i, j] = 0	# Will be visible
    
    if noise:
        sigma = scipy.ndimage.filters.generic_filter(np.ma.filled(data, np.nan), rms_filter, size=(params.stats_bp_window_t, params.stats_bp_window_f))
        sigma = np.ma.filled(sigma, 0)
        print "Creating", os.path.basename(filename)[:-3]+"_"+ant+"_rms.fits"
        f = open(os.path.basename(filename)[:-3]+"_"+ant+"_rms.fits", "wb")
        f.write(fits_header)
        sigma = sigma[::-1]		# have to flip
        sigma.astype(np.float32).byteswap().tofile(f)
        f.close()
        
    if not flag:
        data, bottom, top = scipy.stats.sigmaclip(data, low=params.ds9_clip, high=params.ds9_clip)
        print "Clipping for plotting. Clipping is necessary when flagging is not used,\ndue to extreme peaks, otherwise no detail is visible."
        higher = np.full((T_ant[ant].shape[0], T_ant[ant].shape[1]), top)
        lower = np.full((T_ant[ant].shape[0], T_ant[ant].shape[1]), bottom)
        data = np.where(T_ant[ant]>top, higher, T_ant[ant])
        data = np.where(data<bottom, lower, data)
        
    data = pad_data(data)		# pad out to 2400 channels

    data = np.ma.filled(data, 0)
    
    print "Creating", os.path.basename(filename)[:-3]+"_"+ant+".fits"
    f = open(os.path.basename(filename)[:-3]+"_"+ant+".fits", "wb")
    f.write(fits_header)
    data = data[::-1]		# have to flip
    data.astype(np.float32).byteswap().tofile(f)
    f.close()
    
    if not no_show:
        print "Running DS9"
        os.system("ds9 "+os.path.basename(filename)[:-3]+"_"+ant+".fits")


if __name__ == "__main__":
    import optparse, sys

    usage = '%prog [opts] filename_of_hdf5_observation antenna      # Antenna is like 254B'
    o = optparse.OptionParser()
    o.set_usage(usage)
    o.set_description(__doc__)
    o.add_option('--flag', dest='flag', action='store_true', default=False,
                 help='Apply flagging. Default: False')
    o.add_option('--all_lsts', dest='all_lsts', action='store_true', default=False,
                 help='Include all LSTs, not just when Galaxy and Sun are down. Default: False.')
    o.add_option('--median', dest='median', action='store_true', default=False,
                 help='Use a median filter in the sum threhold flagging. Default: False.')
    o.add_option('--rms', dest='noise', action='store_true', default=False,
                 help='Create a FITs file containing the RMS over the data (after flattening data). The value of a pixel is the RMS at that point. Default: False.')   
    o.add_option('--new_cal', dest='new_cal', action='store_true', default=False,
                 help='Use the 2018 calibration files and method. Default: False.')   

    o.add_option('--no_show', dest='no_show', action='store_true', default=False,
                 help="Don't display the plot on screen using DS9. A FITS file is always created with the same name as the h5 file, but with antenna appended. This can be loaded into another FITS viewer if you don't have DS9. Default: False.")

    opts, args = o.parse_args(sys.argv[1:])

    try:
        filename = args[0]
        ant = args[1]
    except:
        o.print_help()
        exit()

    if opts.new_cal and ant == "252B":
      print "No 2018 calibration data for 252B"
      exit(1)
      
        
    params.median = opts.median
    quicklook(filename, ant, opts.flag,opts.noise, opts.no_show, opts.all_lsts, opts.new_cal)
