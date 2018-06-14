#!/usr/bin/env python

import os
import sys
import ephem
import numpy as np
import tables as tb
from leda_cal.leda_cal import *
from leda_cal.params import params

ovro_location = ('37.2397808', '-118.2816819', 1183.4839)
ovro = ephem.Observer(); (ovro.lat, ovro.lon, ovro.elev) = ovro_location
sun = ephem.Sun()

gal_center = ephem.FixedBody()  
gal_center._ra  =  '17 45 40.04'
gal_center._dec = '-29 00 28.1'
gal_center.name = "Galactic Center"

def main(args, all_lsts):
    filenames = args
    results = {}
    for filename in filenames:
        print "Working on '%s'" % os.path.basename(filename)
        h5 = tb.open_file(filename)

        T_ant = apply_calibration(h5)
        lst_stamps = T_ant['lst']
        if len(lst_stamps) == 0:
          print "No LSTS in file"
          exit(1)

        # Report discontinuities in time
        for i in range(1,len(lst_stamps)):
          if lst_stamps[i]-lst_stamps[i-1] > 1/60.0:	# 1 minute
            print "Discontinuity at LST", lst_stamps[i], (lst_stamps[i]-lst_stamps[i-1])*60*60, "seconds"

        utc_stamps = T_ant['utc']
        # Work out altitude of Gal center and Sun. Use whichever is highest
        # and put that in the padding, which is the stripe.
        unusable_lsts = []
        for i, d in enumerate(utc_stamps):
          ovro.date = d
          sun.compute(ovro)
          gal_center.compute(ovro)
          if sun.alt > params.sun_down*np.pi/180 or gal_center.alt > params.galaxy_down*np.pi/180:
            unusable_lsts.append(i)

        # Delete sun up LSTS
        if not all_lsts:
           print "Cutting out times when sun/galaxy up"
           lst_stamps = np.delete(lst_stamps, unusable_lsts, axis=0)
           utc_stamps = np.delete(utc_stamps, unusable_lsts, axis=0)
           if len(lst_stamps) == 0:
             print "No LSTs available at night time (use --all_lsts to see all)"
             continue
           ylims = ( lst_stamps[0], lst_stamps[-1] )
           print len(lst_stamps), "usable LSTs"
        else: print "Using all LSTs"
        if len(lst_stamps) == 0:
          print "There is no data to display (number of LSTs is 0)"
          continue
        results[os.path.basename(filename)] = lst_stamps

        h5.close()

    filenames = list(results.keys())
    filenames.sort()
    ns = max([len(os.path.splitext(f)[0]) for f in filenames])
    print ''
    for filename in filenames:
        lsts = results[filename]
        tInt = round((lsts[1] - lsts[0])*3600.0)
        nValid = 0.90 * (0.5*3600 / tInt)	# We want the half hour bins to be at least 50% full
        output = ''
        for l in np.linspace(0, 23.5, 48):
            valid = np.where( (lsts>=l) & (lsts<(l+0.5)) )[0]
            if len(valid) == 0:
                output += ' '
            elif len(valid) >= nValid:
                output += '|'
            else:
                output += '/' if output[-1] == ' ' else '\\'
        fmt = "%%%is | 0%%s24" % (ns)
        print fmt % (os.path.splitext(filename)[0], output)

if __name__ == "__main__":
    import optparse, sys

    usage = '%prog [opts] filename_of_hdf5_observation'
    o = optparse.OptionParser()
    o.set_usage(usage)
    o.set_description(__doc__)
    o.add_option('--all_lsts', dest='all_lsts', action='store_true', default=False,
      help='Include all LSTs, not just when Galaxy and Sun are down. A day/night stripe is printed on the right. Default: False.')
    
    opts, args = o.parse_args(sys.argv[1:])
    
    if len(args) < 1:
      o.print_help()
      exit(1)
    else: filenames = args

    main(filenames, opts.all_lsts)
