# LEDA NEC4 Simulations

## Overview
This directory contains a collection of NEC4 simulations for a LWA dipole with a LEDA front end
with a 3 m by 3 m ground screen over a ground with a conductivity of 0.005 S/m

## Files
leda_[xy]e[tp]_1.nec - Base NEC input model for which all other frequencies are derived.  The four files cover 
both polarizations and E and H plane waves.

leda_nec_input.tar.gz - Tarball containing the NEC input models at a variety of frequencies.

leda_nec_output.tar.gz - Tarball containing the NEC output models at a variety of frequencies.

leda_emp_fit_[0-9]+.npz - Python .npz file containing a parameterization following LWA Memo #175 for the NEC 
output model for the corresponding frequency.

leda-dipole-emp.npz - Python .npz file containing the frequency polynomials for the empirical dipole model.

ime.nec - NEC input model for the impedance mis-match efficiency.

ime.out - NEC output model for the impedance mis-match efficiency.

ime.txt - Text version of the dipole IME.

## Scripts
buildModels.py - Script that takes the base NEC models and creates runs at a variety of frequencies.

fitEmp.py - Takes a NEC output file for the dipole power pattern and fits the empirical model from LWA Memo #175.

fitSph.py - Takes a NEC output file for the dipole power pattern and fits spherical harmonics.

updateFits.py - Script that automates running fitEmp.py and fitSph.py.

calcFreq.py - Take the leda_emp_fit_[0-9]+.npz files and fits polynomials and frequencies to each term.

compare74.py - Script that compares the empirical fit with the NEC output model and LWA Memo #175.

extractIME.py - Script that takes a NEC IME output model and converts it to a text file.

imeComparison.py - Script to compare the IME from a NEC output model with the BurnsZ.txt file that contains the IME 
for an LWA dipole as modeled by NRL.
                        



