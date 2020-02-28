#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy

from scipy.optimize import fmin, fmin_cg, fmin_powell, leastsq

from matplotlib import pyplot as plt

def parseNEC(filename):
    theta = []
    phi = []
    amp = []

    fh = open(filename, 'r')
    while True:
        line = fh.readline()
        if line.find('EXCITATION') != -1:
            junk = fh.readline()
            toUse = fh.readline()
            fields = toUse.split()
            theta.append(float(fields[3]))
            phi.append(float(fields[6]))
            while junk.find('AMPS') == -1:
                junk = fh.readline()
            junk = fh.readline()
            toUse = fh.readline()
            fields = toUse.split()
            amp.append(float(fields[6])+1j*float(fields[7]))
        if len(line) == 0:
            break
            
    return numpy.array(theta), numpy.array(phi), numpy.array(amp)

def dipoleModel(B, x):
    za = x[1,:] * numpy.pi/180.0
    az = x[0,:] * numpy.pi/180.0
    
    alpha = B[0]
    beta = B[1]
    if len(B) < 3:
        gamma = 0.0
        delta = 0.0
    elif len(B) < 4:
        delta = 0.0
    else:
        gamma = B[2]
        delta = B[3]
        
    return (1-(2*za/numpy.pi)**alpha)*numpy.cos(za)**beta + gamma*(2*za/numpy.pi)*numpy.cos(za)**delta


def dipoleModelFMin(B, x, y):
    #B = numpy.round(B, 1)
    yFit = dipoleModel(B, x)
    return ((yFit-y)**2).sum()


def dipoleModelLSQ(B, x, y):
    yFit = dipoleModel(B, x)
    return yFit - y


def main(args):
    freq = int(args[0])
    fileBase = ['leda_xep_', 'leda_xet_', 'leda_yep_', 'leda_yet_']
    necData = {}
    for file in fileBase:
        infile = "%s%i.out" % (file, freq)
        junk, pol, junk = infile.split('_', 3)
        theta, phi, amp = parseNEC(infile)
        necData[pol] = {'theta': theta, 'phi': phi, 'amp': amp}
        
    print "X"
    dp = necData['xep']
    dt = necData['xet']
    azaltX = numpy.zeros((2,dp['theta'].size))
    azaltX[0,:] = dp['phi']
    azaltX[1,:] = dp['theta']
    powerX = (numpy.abs(dp['amp'])**2 + numpy.abs(dt['amp'])**2) / 2.0
    scaleX = powerX.max()
    powerX /= scaleX
    
    valid1 = numpy.where( (azaltX[0,:]%360 == 0) | (azaltX[0,:] == 180) )[0]
    valid2 = numpy.where( (azaltX[0,:] == 90) | (azaltX[0,:] == 270) )[0]
    
    outX1 = fmin(dipoleModelFMin, [1,3,1,2], args=(azaltX[:,valid1], powerX[valid1]), maxiter=10000, maxfun=10000, xtol=1e-3, ftol=1e-3)
    #import code
    #code.interact(local=locals())
    mdlX1 = dipoleModel(outX1, azaltX[:,valid1])
    print outX1
    
    to, tc, ti, tm, te = leastsq(dipoleModelLSQ, outX1, args=(azaltX[:,valid1], powerX[valid1]), maxfev=10000, full_output=True)
    tc *= (dipoleModelLSQ(to, azaltX[:,valid1], powerX[valid1])**2).sum() / (len(powerX[valid1])-len(to))
    print 'minimum function value', dipoleModelFMin(to, azaltX[:,valid1], powerX[valid1])
    for i in range(len(to)):
        print 'LSQ', i, to[i], numpy.sqrt(tc[i,i])
        
    print ((mdlX1-powerX[valid1])**2).sum()
    print ((dipoleModel([1,4,1,1.4], azaltX[:,valid1])-powerX[valid1])**2).sum()
    
    outX2 = fmin(dipoleModelFMin, [9,1], args=(azaltX[:,valid2], powerX[valid2]), maxiter=10000, maxfun=10000, xtol=1e-3, ftol=1e-3)
    mdlX2 = dipoleModel(outX2, azaltX[:,valid2])
    mdlX2 /= mdlX2.max()
    print outX2
    
    ti = [outX2[0], outX2[1], 0.0, 0.0]
    to, tc, ti, tm, te = leastsq(dipoleModelLSQ, ti, args=(azaltX[:,valid2], powerX[valid2]), maxfev=10000, full_output=True)
    tc *= (dipoleModelLSQ(to, azaltX[:,valid2], powerX[valid2])**2).sum() / (len(powerX[valid2])-len(to))
    print 'minimum function value', dipoleModelFMin(to, azaltX[:,valid2], powerX[valid2])
    for i in range(len(to)):
        print 'LSQ', i, to[i], numpy.sqrt(tc[i,i])
        
    print ((mdlX2-powerX[valid2])**2).sum()
    print ((dipoleModel([9,1.2,0,0], azaltX[:,valid2])-powerX[valid2])**2).sum()
    
    
    print "Y"
    dp = necData['yep']
    dt = necData['yet']
    azaltY = numpy.zeros((2,dp['theta'].size))
    azaltY[0,:] = dp['phi']
    azaltY[1,:] = dp['theta']
    powerY = (numpy.abs(dp['amp'])**2 + numpy.abs(dt['amp'])**2) / 2.0
    scaleY = powerY.max()
    powerY /= scaleY
    
    valid3 = numpy.where( (azaltY[0,:]%360 == 0) | (azaltY[0,:] == 180) )[0]
    valid4 = numpy.where( (azaltY[0,:] == 90) | (azaltY[0,:] == 270) )[0]
    
    outY1 = fmin(dipoleModelFMin, [9,1], args=(azaltY[:,valid3], powerY[valid3]), maxiter=10000, maxfun=10000, xtol=1e-3, ftol=1e-3)
    mdlY1 = dipoleModel(outY1, azaltY[:,valid3])
    mdlY1 /= mdlY1.max()
    print outY1
    
    print ((mdlY1-powerY[valid3])**2).sum()
    print ((dipoleModel([9,1.2,0,0], azaltY[:,valid3])-powerY[valid3])**2).sum()
    
    outY2 = fmin(dipoleModelFMin, [1,3,1,2], args=(azaltY[:,valid4], powerY[valid4]), maxiter=10000, maxfun=10000, xtol=1e-3, ftol=1e-3)
    mdlY2 = dipoleModel(outY2, azaltY[:,valid4])
    mdlY2 /= mdlY2.max()
    print outY2
    
    print ((mdlY2-powerY[valid4])**2).sum()
    print ((dipoleModel([1,4,1,1.4], azaltY[:,valid4])-powerY[valid4])**2).sum()
    
    numpy.savez('leda_emp_fit_%i.npz' % freq, scaleX=scaleX, scaleY=scaleY, out1=outX1, out2=outX2, out3=outY1, out4=outY2)


if __name__ == "__main__":
    numpy.seterr(all='ignore')
    main(sys.argv[1:])
    
