#!/usr/bin/env python

from __future__ import print_function

import os
import re
import sys
import glob
import numpy
import argparse


_FILENAME_RE = re.compile('driftcurve_(?P<site>[a-z0-9]+)_(?P<pol>(EW|NS))_(?P<freq>\d+(\.\d*)?)_(?P<model>[a-z]+)_(?P<ant>[a-z]+)\.txt')


def main(args):
    globex = 'driftcurve_ovro_%s_*_gsm_nec.txt' % args.polarization.upper()
    if args.lfsm:
        globex = globex.replace('_gsm', '_lfsm')
    if args.empirical:
        globex = globex.replace('_nec', '_emp')
    filenames = glob.glob(globex)
    filenames.sort()
    
    freqs = []
    data = []
    for filename in filenames:
        # Make sure the frequence is in range before we load the file
        mtch = _FILENAME_RE.match(filename)
        freq = float(mtch.group('freq'))
        if freq < args.freq_min or freq > args.freq_max:
            continue
            
        # Load the data and save it
        subdata = numpy.loadtxt(filename)
        freqs.append(freq)
        data.append(subdata)
    freqs = numpy.array(freqs)
    data = numpy.array(data)
    # Add in the temperature offset
    data[:,:,1] += args.offset
    
    print("# Sky Model: %s" % ('LFSM' if args.lfsm else 'GSM',))
    print("# Antenna Model: %s" % ('empirical' if args.empirical else 'NEC',))
    print("# Polarization: %s" % args.polarization)
    print("# Temperature Offset: %.3f K" % args.offset)
    print("# Fit Frequency Range: %.3f to %.3f MHz" % (args.freq_min, args.freq_max))
    #print("#                      (Used %s MHz)" % ','.join(["%.3f" % f for f in freqs]))
    print("# Reference Frequency: %.3f MHz" % args.freq_ref)
    print("#")
    print("# Columns:  LST, spectra index, T0")
    for i in range(data.shape[1]):
        lst = data[0,i,0]
        subdata = data[:,i,1]
        fit = numpy.polyfit(numpy.log10(freqs/args.freq_ref), numpy.log10(subdata), 1)
        print("%.6f  %.6f  %.3f" % (lst, fit[0], 10**fit[1]))


if __name__ == "__main__":
    def polarization(value):
        if value.upper() not in ('EW', 'NS'):
            raise argparse.ArgumentError
        return value.upper()
    
    parser = argparse.ArgumentParser(
        description='calculate the spectral index of the model sky as a function of LST and write the resutls to stdout',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('-p', '--polarization', type=polarization, default='EW',
                        help='polarization of the observations (NS or EW)')
    parser.add_argument('-e', '--empirical', action='store_true',
                        help='enable empirical corrections to the dipole model (valid from 35 to 80 MHz)')
    parser.add_argument('-l', '--lfsm', action='store_true',
                        help='use LFSM instead of GSM')
    parser.add_argument('-o', '--offset', type=float, default=0.0,
                        help='temperature offset to apply before computing')
    parser.add_argument('-n', '--freq-min', type=float, default=40.0,
                        help='minimum frequency in MHz to use for fitting')
    parser.add_argument('-x', '--freq-max', type=float, default=80.0,
                        help='maximum frequency in MHz to use for fitting')
    parser.add_argument('-r', '--freq-ref', type=float, default=60.0,
                        help='reference frequency in MHz of the spectra index fit')
    args = parser.parse_args()
    main(args)
    