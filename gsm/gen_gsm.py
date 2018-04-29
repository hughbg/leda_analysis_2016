"""
gen_gsm.py
----------

Convert driftscan.py outputs into a single hkl file.

Output file is 15 x 144 in size, one file per polarization.
First col is lst, all others are frequencies 30 - 90 MHz.

"""

import numpy as np
import hickle as hkl

files_ew = [
    'driftcurve_lwa1_EW_30.00.txt',
    'driftcurve_lwa1_EW_35.00.txt',
    'driftcurve_lwa1_EW_40.00.txt',
    'driftcurve_lwa1_EW_45.00.txt',
    'driftcurve_lwa1_EW_50.00.txt',
    'driftcurve_lwa1_EW_55.00.txt',
    'driftcurve_lwa1_EW_60.00.txt',
    'driftcurve_lwa1_EW_65.00.txt',
    'driftcurve_lwa1_EW_70.00.txt',
    'driftcurve_lwa1_EW_75.00.txt',
    'driftcurve_lwa1_EW_80.00.txt',
    'driftcurve_lwa1_EW_85.00.txt',
    'driftcurve_lwa1_EW_90.00.txt'
]

files_ns = [
    'driftcurve_lwa1_NS_30.00.txt',
    'driftcurve_lwa1_NS_35.00.txt',
    'driftcurve_lwa1_NS_40.00.txt',
    'driftcurve_lwa1_NS_45.00.txt',
    'driftcurve_lwa1_NS_50.00.txt',
    'driftcurve_lwa1_NS_55.00.txt',
    'driftcurve_lwa1_NS_60.00.txt',
    'driftcurve_lwa1_NS_65.00.txt',
    'driftcurve_lwa1_NS_70.00.txt',
    'driftcurve_lwa1_NS_75.00.txt',
    'driftcurve_lwa1_NS_80.00.txt',
    'driftcurve_lwa1_NS_85.00.txt',
    'driftcurve_lwa1_NS_90.00.txt'
]


gsm_ew = np.zeros((144, 14))
gsm_ns = np.zeros((144, 14))

ii = 0
for filename in files_ew:
    a = np.genfromtxt(filename)
    if ii == 0:
        gsm_ew[:, 0] = a[:, 0]
        gsm_ew[:, 1] = a[:, 1]
        ii = 2
    else:
        gsm_ew[:, ii] = a[:, 1]
        ii += 1

ii = 0
for filename in files_ns:
    a = np.genfromtxt(filename)
    if ii == 0:
        gsm_ns[:, 0] = a[:, 0]
        gsm_ns[:, 1] = a[:, 1]
        ii = 2
    else:
        gsm_ns[:, ii] = a[:, 1]
        ii += 1


hkl.dump(gsm_ew, 'gsm_ew.hkl')
hkl.dump(gsm_ns, 'gsm_ns.hkl')

