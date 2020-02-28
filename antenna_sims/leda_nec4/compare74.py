#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import sys
import numpy

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


def steve74E(theta):
	alpha = 1.0
	beta  = 4.0
	gamma = 1.0
	delta = 1.4
	
	out = (1-(2*theta/numpy.pi)**alpha)*numpy.cos(theta)**beta
	out += gamma*(2*theta/numpy.pi)*numpy.cos(theta)**delta
	return out


def steve74H(theta):
	alpha = 9.0
	beta  = 1.2
	gamma = 0.0
	delta = 0.0
	
	out = (1-(2*theta/numpy.pi)**alpha)*numpy.cos(theta)**beta
	out += gamma*(2*theta/numpy.pi)*numpy.cos(theta)**delta
	return out


def model(params, theta):
	alpha = params[0]
	beta  = params[1]
	if len(params) < 3:
		gamma = 0.0
		delta = 0.0
	elif len(params) < 4:
		gammas = params[2]
		delta = 0.0
	else:
		gamma = params[2]
		delta = params[3]
	
	out = (1-(2*theta/numpy.pi)**alpha)*numpy.cos(theta)**beta
	out += gamma*(2*theta/numpy.pi)*numpy.cos(theta)**delta
	return out


def main(args):
	freq = 74
	fileBase = ['leda_xep_', 'leda_xet_', 'leda_yep_', 'leda_yet_']
	necData = {}
	for file in fileBase:
		infile = "%s%i.out" % (file, freq)
		junk, pol, junk = infile.split('_', 3)
		theta, phi, amp = parseNEC(infile)
		necData[pol] = {'theta': theta, 'phi': phi, 'amp': amp}

	dp = necData['xep']
	dt = necData['xet']
	azaltX = numpy.zeros((2,dp['theta'].size))
	azaltX[0,:] = dp['phi']
	azaltX[1,:] = dp['theta']
	powerX = (numpy.abs(dp['amp'])**2 + numpy.abs(dt['amp'])**2) / 2.0
	powerX /= powerX.max()

	valid1 = numpy.where( (azaltX[0,:]%360 == 0) | (azaltX[0,:] == 180) )[0]
	valid2 = numpy.where( (azaltX[0,:] == 90) | (azaltX[0,:] == 270) )[0]
	
	dataDict = numpy.load('leda_emp_fit_74.npz')
	fitE = dataDict['out1']
	fitH = dataDict['out2']
	
	za = numpy.linspace(0,numpy.pi/2, 181)
	fig = plt.figure()
	ax1 = fig.add_subplot(1, 1, 1)
	#ax1 = fig.add_subplot(2, 1, 1)
	#ax2 = fig.add_subplot(2, 1, 2)
	
	ax1.plot(azaltX[1,valid1], numpy.log10(powerX[valid1])*10, color='blue', marker='x', linestyle=' ', label='E-plane')
	ax1.plot(azaltX[1,valid2], numpy.log10(powerX[valid2])*10, color='green', marker='x', linestyle=' ', label='H-plane')
	
	ax1.plot(za*180/numpy.pi, numpy.log10(model(fitE, za))*10, color='blue', label='E-plane LEDA')
	ax1.plot(za*180/numpy.pi, numpy.log10(model(fitH, za))*10, color='green', label='E-plane LEDA')
	
	ax1.plot(za*180/numpy.pi, numpy.log10(steve74E(za))*10, color='blue', linestyle='--', label='E-plane LWA 175')
	ax1.plot(za*180/numpy.pi, numpy.log10(steve74H(za))*10, color='green', linestyle='--', label='H-plane LWA 175')
	
	#ax2.plot(azaltX[1,valid1], 10**(numpy.log10(powerX[valid1])-numpy.log10(model(fitE, azaltX[1,valid1]*numpy.pi/180))), color='blue')
	#ax2.plot(azaltX[1,valid1], 10**(numpy.log10(powerX[valid2])-numpy.log10(model(fitH, azaltX[1,valid1]*numpy.pi/180))), color='green')
	
	#ax2.plot(azaltX[1,valid1], 10**(numpy.log10(powerX[valid1])-numpy.log10(steve74E(azaltX[1,valid1]*numpy.pi/180))), color='blue', linestyle='--')
	#ax2.plot(azaltX[1,valid1], 10**(numpy.log10(powerX[valid2])-numpy.log10(steve74H(azaltX[1,valid1]*numpy.pi/180))), color='green', linestyle='--')
	
	leg = ax1.legend(loc=0)
	for l in leg.get_lines():
		l.set_linewidth(1.2)
	ax1.set_xlabel('Zenith Angle [deg]')
	ax1.set_ylabel('Normalized Power Pattern [dB]')
	ax1.set_ylim([-25, 2])
	#ax2.set_ylim([0.8,1.2])
	ax1.grid()
	plt.show()
	
	fig.savefig('comparison-74MHz.png')


if __name__ == "__main__":
	main(sys.argv[1:])
