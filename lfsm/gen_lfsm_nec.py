"""
gen_lfsm.py
----------

Run driftscan.py and convert the outputs into a single hkl file.

Output file is 15 x 144 in size, one file per polarization.
First col is lst, all others are frequencies 30 - 90 MHz.

"""

import os
import numpy as np
import hickle as hkl

for freq in xrange(30, 95, 5):
    os.system("python driftcurve.py -o -l -f %i -p EW" % freq)
    os.system("python driftcurve.py -o -l -f %i -p NS" % freq)

files_ew = [
    'driftcurve_ovro_EW_30.00.txt',
    'driftcurve_ovro_EW_35.00.txt',
    'driftcurve_ovro_EW_40.00.txt',
    'driftcurve_ovro_EW_45.00.txt',
    'driftcurve_ovro_EW_50.00.txt',
    'driftcurve_ovro_EW_55.00.txt',
    'driftcurve_ovro_EW_60.00.txt',
    'driftcurve_ovro_EW_65.00.txt',
    'driftcurve_ovro_EW_70.00.txt',
    'driftcurve_ovro_EW_75.00.txt',
    'driftcurve_ovro_EW_80.00.txt',
    'driftcurve_ovro_EW_85.00.txt',
    'driftcurve_ovro_EW_90.00.txt'
]

files_ns = [
    'driftcurve_ovro_NS_30.00.txt',
    'driftcurve_ovro_NS_35.00.txt',
    'driftcurve_ovro_NS_40.00.txt',
    'driftcurve_ovro_NS_45.00.txt',
    'driftcurve_ovro_NS_50.00.txt',
    'driftcurve_ovro_NS_55.00.txt',
    'driftcurve_ovro_NS_60.00.txt',
    'driftcurve_ovro_NS_65.00.txt',
    'driftcurve_ovro_NS_70.00.txt',
    'driftcurve_ovro_NS_75.00.txt',
    'driftcurve_ovro_NS_80.00.txt',
    'driftcurve_ovro_NS_85.00.txt',
    'driftcurve_ovro_NS_90.00.txt'
]


lfsm_ew = np.zeros((144, 14))
lfsm_ns = np.zeros((144, 14))

ii = 0
for filename in files_ew:
    a = np.genfromtxt(filename)
    if ii == 0:
        lfsm_ew[:, 0] = a[:, 0]
        lfsm_ew[:, 1] = a[:, 1]
        ii = 2
    else:
        lfsm_ew[:, ii] = a[:, 1]
        ii += 1

ii = 0
for filename in files_ns:
    a = np.genfromtxt(filename)
    if ii == 0:
        lfsm_ns[:, 0] = a[:, 0]
        lfsm_ns[:, 1] = a[:, 1]
        ii = 2
    else:
        lfsm_ns[:, ii] = a[:, 1]
        ii += 1


hkl.dump(lfsm_ew, 'lfsm_ew.hkl')
hkl.dump(lfsm_ns, 'lfsm_ns.hkl')

