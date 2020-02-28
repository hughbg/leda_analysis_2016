#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import numpy

from scipy.interpolate import interp1d

from lsl.sim.nec_util import NECImpedance

filename = sys.argv[1]
imp = NECImpedance(filename)
print("""# Columns:
#  1) Frequency [MHz]
#  2) Impedance - Real [ohm]
#  3) Impedance - Imaginary [ohm]
#  4) Impedance Mis-Match Efficiency (1-|Gamma|^2)""")
for i in xrange(imp.freqs.size):
    z = imp.z[i]
    fe = 200.0      # 200 Ohm for a LEDA front end
    g = (z - fe) / (z + fe)
    ime = 1 - numpy.abs(g)**2
    print("%-3i  %8.3f  %8.3f  %11.9f" % (imp.freqs[i], z.real, z.imag, ime))
