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
from leda_cal.useful import statistics

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

def quicklook(filename, ant, flag, no_show, all_lsts):
    global fits_header

    h5 = tb.open_file(filename)
    T_ant = apply_calibration(h5)
    f_leda = T_ant['f']
    lst_stamps = T_ant['lst']
    utc_stamps = T_ant['utc']

    # Report discontinuities in time
    for i in range(1,len(lst_stamps)):
      if lst_stamps[i]-lst_stamps[i-1] > 1/60.0:	# 1 minute
        print "Discontinuity at LST", lst_stamps[i], (lst_stamps[i]-lst_stamps[i-1])*60*60, "seconds"

    # Work out altitude of Gal center and Sun. Use whichever is highest
    # and put that in the padding, which is the stripe.
    unusable_lsts = []
    for i, d in enumerate(utc_stamps):
      ovro.date = d
      sun.compute(ovro)
      gal_center.compute(ovro)
      if sun.alt > -15*np.pi/180 or gal_center.alt > -15*np.pi/180:
        unusable_lsts.append(i)
  
    # Delete sun up LSTS
    if not all_lsts:
       print "Cutting out Sun/Galaxy up. LSTS in DS9 will be WRONG!!!!!!"
       lst_stamps = np.delete(lst_stamps, unusable_lsts, axis=0)
       utc_stamps = np.delete(utc_stamps, unusable_lsts, axis=0)
       print len(lst_stamps), "usable LSTs"
    else: print "Using all LSTs"
    if len(lst_stamps) == 0:
      print "There is no data to display (number of LSTs is 0)"
      exit(1) 
    
    fits_header = fits_header.replace("XXXX", ( "%04d" % T_ant[ant].shape[1] ))
    fits_header = fits_header.replace("YYYY", ( "%04d" % len(lst_stamps) ))

    lst_diff = str((lst_stamps[0]-lst_stamps[-1])/(len(lst_stamps)-1))[:10]
    freq_diff = str((f_leda[-1]-f_leda[0])/(len(f_leda)-1))[:10]
    fits_header = fits_header.replace("LLLLLLLLLL", lst_diff)
    fits_header = fits_header.replace("FFFFFFFFFF", freq_diff)

    lst_start = str(lst_stamps[-1])[:6]			# start and the end for FITs in our orientation
    freq_start = str(f_leda[0])[:6]
    fits_header = fits_header.replace("STARTL", lst_start)
    fits_header = fits_header.replace("STARTF", freq_start)
    
   

    if flag:
      print "Flagging"
      data = rfi_flag(T_ant[ant], thr_f=0.2, thr_t=0.2, rho=1.5,
               bp_window_f=16, bp_window_t=16,
               max_frac_f=0.5, max_frac_t=0.5, freqs=f_leda)
    else:   
      data, bottom, top = scipy.stats.sigmaclip(T_ant[ant])
      print "Clipping for plotting"
      higher = np.full((T_ant[ant].shape[0], T_ant[ant].shape[1]), top)
      lower = np.full((T_ant[ant].shape[0], T_ant[ant].shape[1]), bottom)

      data = np.where(T_ant[ant]>top, higher, T_ant[ant])
      data = np.where(data<bottom, lower, data)

    if not all_lsts:  
      if flag: mask = np.delete(data.mask, unusable_lsts, axis=0)
      data = np.delete(data, unusable_lsts, axis=0)
      if flag: data.mask = mask


    statistics(data)

    # Test the image is oriented right
    #for i in range(200):
    #  for j in range(200): data[i, j] = 0	# Will be visible
	
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
    o.add_option('--no_show', dest='no_show', action='store_true', default=False,
      help="Don't display the plot on screen using DS9. A FITS file is always created with the same name as the h5 file, but with antenna appended. This can be loaded into another FITS viewer if you don't have DS9. Default: False.")

    opts, args = o.parse_args(sys.argv[1:])

    try:
        filename = args[0]
        ant = args[1]
    except:
        o.print_help()
        exit()
    
    quicklook(filename, ant, opts.flag, opts.no_show, opts.all_lsts)
  

    