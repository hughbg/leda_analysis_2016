#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import glob
import numpy
from lsl.statistics import robust
from matplotlib import pyplot as plt

def mainSph(args):
	files = glob.glob('leda_sph_fit_*.npz')
	files = sorted(files)
	nFreq = len(files)
	print nFreq
	
	dataDict = numpy.load(files[0])
	nTerms = len(dataDict['termsX'])
	
	freq = numpy.zeros(nFreq)
	x = numpy.zeros((nFreq, nTerms), dtype=numpy.complex64)
	y = numpy.zeros((nFreq, nTerms), dtype=numpy.complex64)
	a = numpy.zeros((nFreq, nTerms), dtype=numpy.complex64)
	
	i = 0
	for file in files:
		dataDict = numpy.load(file)
		f, junk = file.split('.')
		junk, junk, f = f.rsplit('_', 2)
		f = int(f)

		try:
			x[i,:] = dataDict['termsX']
			y[i,:] = dataDict['termsY']
			a[i,:] = dataDict['termsAIPY']
		except:
			pass
		freq[i] = f
		i = i + 1

	tooSmall = numpy.where( numpy.abs(x) < 2e-7 )
	x[tooSmall] = 0
	tooSmall = numpy.where( numpy.abs(y) < 2e-7 )
	y[tooSmall] = 0
	tooSmall = numpy.where( numpy.abs(a) < 2e-7 )
	a[tooSmall] = 0

	fitAipy = numpy.zeros((8, nTerms), dtype=numpy.complex_)
	fitX = numpy.zeros((8, nTerms), dtype=numpy.complex_)
	fitY = numpy.zeros((8, nTerms), dtype=numpy.complex_)
	t = 0
	for t in range(nTerms):
		fitAipy[:,t].real = robust.polyfit(freq/1000.0, a[:,t].real, 7)
		fitAipy[:,t].imag = robust.polyfit(freq/1000.0, a[:,t].imag, 7)
		fitX[:,t].real = robust.polyfit(freq*1e6, x[:,t].real, 7)
		fitX[:,t].imag = robust.polyfit(freq*1e6, x[:,t].imag, 7)
		fitY[:,t].real = robust.polyfit(freq*1e6, y[:,t].real, 7)
		fitY[:,t].imag = robust.polyfit(freq*1e6, y[:,t].imag, 7)
		t += 1

	numpy.savez('beam-shape.npz', coeffs=fitAipy)
	numpy.savez('leda-dipole-sph.npz', fitX=fitX, fitY=fitY)

	t = 0
	l = 0
	m = 0
	for p in xrange(int(numpy.ceil(x.shape[1]/25))+1):
		fig = plt.figure()
		for i in xrange(25):
			if t >= x.shape[1]:
				break

			ax = fig.add_subplot(5, 5, i+1)
			
			ax.plot(freq, x[:,t].real, marker='x', linestyle=' ', color='blue')
			ax.plot(freq, x[:,t].imag, marker='o', linestyle=' ', color='blue')
			ax.plot(freq, y[:,t].real, marker='x', linestyle=' ', color='red')
			ax.plot(freq, y[:,t].imag, marker='o', linestyle=' ', color='red')
			ax.plot(numpy.linspace(10,88,100), numpy.polyval(fitX[:,t].real, numpy.linspace(10,88,100)*1e6), linestyle='-', color='blue')
			ax.plot(numpy.linspace(10,88,100), numpy.polyval(fitX[:,t].imag, numpy.linspace(10,88,100)*1e6), linestyle='-.', color='blue')
			ax.plot(numpy.linspace(10,88,100), numpy.polyval(fitY[:,t].real, numpy.linspace(10,88,100)*1e6), linestyle='-', color='red')
			ax.plot(numpy.linspace(10,88,100), numpy.polyval(fitY[:,t].imag, numpy.linspace(10,88,100)*1e6), linestyle='-.', color='red')

			ax.set_title('l=%i, m=%i' % (l, m))
			m = m + 1
			if m > l:
				l = l + 1
				m = 0
	
			t += 1

	plt.show()


def mainEmp(args):
	files = glob.glob('leda_emp_fit_*.npz')
	files = sorted(files)
	nFreq = len(files)
	print nFreq
	
	freq = numpy.zeros(nFreq)
	x1 = numpy.zeros((nFreq,4))
	x2 = numpy.zeros((nFreq,4))
	y1 = numpy.zeros((nFreq,4))
	y2 = numpy.zeros((nFreq,4))
	
	i = 0
	for file in files:
		dataDict = numpy.load(file)
		f, junk = file.split('.')
		junk, junk, f = f.rsplit('_', 2)
		f = int(f)

		freq[i] = f
		x1[i,:] = dataDict['out1']
		if len(dataDict['out2']) == 4:
			x2[i,:] = dataDict['out2']
		else:
			x2[i,0:2] = dataDict['out2']
		y1[i,:] = dataDict['out4']
		if len(dataDict['out3']) == 4:
			y2[i,:] = dataDict['out3']
		else:
			y2[i,0:2] = dataDict['out3']
		i = i + 1
	
	fitX = numpy.zeros((2, 4, 12))
	fitY = numpy.zeros((2, 4, 12))

	fig = plt.figure()
	label = ['$\\alpha$', '$\\beta$', '$\\gamma$', '$\\delta$']
	axs = []
	for t in range(4):
		ax = fig.add_subplot(2, 2, t+1)
		axs.append(ax)
		#ax.plot(freq, x1[:,t], marker='x', linestyle=' ', color='blue')
		ax.plot(freq, y1[:,t], marker='x', linestyle=' ', color='blue')
		#ax.plot(numpy.linspace(10,88,100), numpy.polyval(robust.polyfit(freq, x1[:,t], 9), numpy.linspace(10,88,100)), color='blue')
		#ax.plot(numpy.linspace(10,88,100), numpy.polyval(robust.polyfit(freq, y1[:,t], 9), numpy.linspace(10,88,100)), color='blue')
		ax.set_title('%s' % label[t])

		fitX[0,t,:] = numpy.polyfit(freq*1e6, x1[:,t], 11)
		fitY[1,t,:] = numpy.polyfit(freq*1e6, y1[:,t], 11)
	#plt.draw()

	#fig = plt.figure()
	for t in range(4):
		ax = axs[t]
		#ax.plot(freq, x2[:,t], marker='x', linestyle=' ', color='green')
		ax.plot(freq, y2[:,t], marker='x', linestyle=' ', color='green')
		#ax.plot(numpy.linspace(10,88,100), numpy.polyval(robust.polyfit(freq, x2[:,t], 9), numpy.linspace(10,88,100)), color='green')
		#ax.plot(numpy.linspace(10,88,100), numpy.polyval(robust.polyfit(freq, y2[:,t], 9), numpy.linspace(10,88,100)), color='green')
		#ax.set_title('Part: 2, Term: %i' % t)

		fitX[1,t,:] = numpy.polyfit(freq*1e6, x2[:,t], 11)
		fitY[0,t,:] = numpy.polyfit(freq*1e6, y2[:,t], 11)
		
	for ax in axs:
		ax.set_xlabel('Freq. [MHz]')
	plt.draw()

	plt.show()
	
	numpy.savez('leda-dipole-emp.npz', fitX=fitX, fitY=fitY)

if __name__ == "__main__":
	mainEmp(sys.argv[1:])
