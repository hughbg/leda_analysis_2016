#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import aipy
import numpy

from lsl.misc import mathutil

from matplotlib import pyplot as plt

def parseNEC2(filename):
    output = numpy.zeros((360,91), dtype=numpy.complex64)
    
    fh = open(filename, 'r')
    
    n = 0
    while True:
        line = fh.readline()
        if line.find('EXCITATION') >= 0:
            parts = []
            for l in xrange(12):
                parts.append( fh.readline() )
            fieldsAngle = parts[1].split()
            fieldsCurrent = parts[11].split()

            theta = 90.0 - float(fieldsAngle[3])
            phi = float(fieldsAngle[6])
            if theta < 0 or phi > 359:
                continue
                
            output[int(phi), int(theta)] = float(fieldsCurrent[6]) + 1j*float(fieldsCurrent[7])
            
        if len(line) == 0:
            break
    return output


def fit(az, alt, data, lmax, label=None, degrees=True, realOnly=True):
    """Fit a given set of data using spherical harmonics.  Return the list
    of coefficients."""
    
    terms = mathutil.sphfit(az, alt, data, lmax=lmax, degrees=degrees, realOnly=realOnly)
    fit = mathutil.sphval(terms, az, alt, degrees=degrees, realOnly=realOnly)
    diff = data - fit
    
    if label is not None:
        print " "+str(label)
    print "  Peak Differences:", data.max(), fit.max()
    print "  Model Differences:", diff.min(), diff.mean(), diff.max()
    print "  Model RMS:", (diff**2).sum()
    
    return terms


def main(args):
    freq = int(args[0])
    print "Working on %.2f MHz" % freq
    
    extP = parseNEC2('leda_xep_%i.out' % freq)
    extT = parseNEC2('leda_xet_%i.out' % freq)
    extX = (numpy.abs(extP)**2 + numpy.abs(extT)**2) / 2.0
    sclX = extX.max()
    extX /= sclX
    
    extP = parseNEC2('leda_yep_%i.out' % freq)
    extT = parseNEC2('leda_yet_%i.out' % freq)
    extY = (numpy.abs(extP)**2 + numpy.abs(extT)**2) / 2.0
    sclY = extY.max()
    extY /= sclY
    
    extAz = numpy.zeros((360,91))
    extAlt = numpy.zeros((360,91))
    for i in xrange(360):
        extAz[i,:] = i
    for i in xrange(91):
        extAlt[:,i] = i
        
    top = aipy.coord.azalt2top(numpy.array([[extAz*numpy.pi/180], [extAlt*numpy.pi/180]]))
    theta, phi = aipy.coord.xyz2thphi(top)
    theta = theta.squeeze()
    phi = phi.squeeze()
    
    l = 12
    from multiprocessing import Pool
    taskPool = Pool(3)
    taskList = []
    
    task = taskPool.apply_async(fit, args=(extAz, extAlt, extX, l), kwds={'label': "X", 'degrees': True, 'realOnly':True})
    taskList.append(('X', task))
    task = taskPool.apply_async(fit, args=(extAz, extAlt, extY, l), kwds={'label': "Y", 'degrees': True, 'realOnly':True})
    taskList.append(('Y', task))
    task = taskPool.apply_async(fit, args=(phi, theta-numpy.pi/2, extX, l), kwds={'label': "AIPY", 'degrees': False, 'realOnly':True})
    taskList.append(('A', task))
    
    taskPool.close()
    taskPool.join()
    
    for code,task in taskList:
        if code == "X":
            termsX = task.get()
        elif code == "Y":
            termsY = task.get()
        else:
            termsA = task.get()
            
    numpy.savez('leda_sph_fit_%i.npz' % freq, l=l, realOnly=True, scaleX=sclX, scaleY=sclY, termsX=termsX, termsY=termsY, termsAIPY=termsA)


if __name__ == "__main__":
    numpy.seterr(all='ignore')
    main(sys.argv[1:])
    
