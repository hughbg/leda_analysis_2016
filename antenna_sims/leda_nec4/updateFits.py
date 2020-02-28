#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Walk through the avaliable NEC files in the current directory and fit those
that appear to have all four of the NEC output files.
"""

import os
import re
import sys
import glob
import getopt


def usage(exitCode=None):
    print """updateFits.py - Compute fits for the various NEC models that have
been computed.

-f, --force               Force compting all fits even if the NPZ file exists
-e, --empirical           Compute only the empirical fits           
-s, --spherical           Compute only the spherical harmonic fits
"""
    
    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseOptions(args):
    config = {}
    # Command line flags - default values
    config['force'] = False
    config['emp'] = True
    config['sph'] = True
    
    # Read in and process the command line flags
    try:
        opts, args = getopt.getopt(args, "hfes", ["help", "force", "empirical", "spherical"])
    except getopt.GetoptError, err:
        # Print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage(exitCode=2)
    
    # Work through opts
    for opt, value in opts:
        if opt in ('-h', '--help'):
            usage(exitCode=0)
        elif opt in  ('-f', '--force'):
            config['force'] = True
        elif opt in ('-e', '--empirical'):
            config['sph'] = False
        elif opt in ('-s', '--spherical'):
            config['emp'] = False
        else:
            assert False
            
    # Return configuration
    return config


def main(args):
    config = parseOptions(args)
    
    necFiles = glob.glob('leda_xep_*.nec')
    necFiles.sort()
    
    readyFreqs = []
    necRE = re.compile(r'leda_(?P<type>[xy]e[tp])_(?P<freq>\d+).nec')
    for necFile in necFiles:
        mtch = necRE.match(necFile)
        freq = int(mtch.group('freq'))
        
        outFiles = 0
        for sType in ['xet', 'xep', 'yet', 'yep']:
            if os.path.exists('leda_%s_%i.out' % (sType, freq)):
                outFiles += 1
                
        if outFiles == 4:
            readyFreqs.append(freq)
            
    nEmp = 0
    nSph = 0
    for freq in readyFreqs:
        if config['emp']:
            if os.path.exists('leda_emp_fit_%i.npz' % freq):
                if config['force']:
                    os.unlink('leda_emp_fit_%i.npz' % freq)
                    exitStatus = os.system('fitEmp.py %i' % freq)
                    nEmp += 1
                else:
                    continue
            else:
                exitStatus = os.system('fitEmp.py %i' % freq)
                nEmp += 1
                
        if config['sph']:
            if os.path.exists('leda_sph_fit_%i.npz' % freq):
                if config['force']:
                    os.unlink('leda_sph_fit_%i.npz' % freq)
                    exitStatus = os.system('fitSph.py %i' % freq)
                    nSph += 1
                else:
                    continue
            else:
                exitStatus = os.system('fitSph.py %i' % freq)
                nSph += 1
                
    print "Ran %i empirical fits and %i spherical fits" % (nEmp, nSph)


if __name__ == "__main__":
    main(sys.argv[1:])
    
