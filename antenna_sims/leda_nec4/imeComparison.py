#!/usr/bin/env python

import os
import sys
import numpy

from scipy.interpolate import interp1d

from lsl.common.paths import data
from lsl.sim.nec_util import NECImpedance


ime0 = numpy.loadtxt(os.path.join(data, 'BurnsZ.txt'))
freq0, ime0 = ime0[:,0], ime0[:,3]

import pylab
pylab.plot(freq0, ime0, linestyle='-', marker='+', label='BurnsZ')
for filename in sys.argv[1:]:
    try:
        imp = NECImpedance(filename)
    except (RuntimeError, ValueError):
        continue
        
    ime = []
    for i in xrange(imp.freqs.size):
        z = imp.z[i]
        fe = 200.0      # 200 Ohm for a LEDA front end
        g = (z - fe) / (z + fe)
        ime.append(1 - numpy.abs(g)**2)
        
    pylab.plot(imp.freqs, ime, label="%s" % os.path.basename(filename))
pylab.legend(loc=0)
pylab.xlabel('Frequency [MHz]')
pylab.ylabel('1-|$\\Gamma$|$^2$')
pylab.show()


