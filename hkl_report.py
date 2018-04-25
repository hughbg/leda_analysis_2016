#!/usr/bin/env python

import hickle, sys, numpy as np
import matplotlib.pyplot as plt

def load_hkl(filename):
  return hickle.load(filename)

def report_hkl(data):
  print "Items present:"
  for key in sorted(data.keys()): 
    try:
      if len(data[key].shape) == 1: print "Key", '"'+key+'"', "Type", type(data[key]), "Element Type", type(data[key][0]), "Shape", data[key].shape
      else: print "Key", '"'+key+'"', "Type", type(data[key]), "Element Type", type(data[key][0][0]), "Shape", data[key].shape, "(time, chan)"
    except: pass
  print "\nStart time", data["utcs"][0]
  print "Options", data["options"]

  # Sometimes you have to replace masked data with a value
  #plt.imshow(np.ma.filled(data["252A"], 0), clim=(1000,10000), aspect="auto")
  #plt.show()


if len(sys.argv) == 1:
  print "Need 1 hkl file as argument"
  exit(1)

d = load_hkl(sys.argv[1])
report_hkl(d)

# Example code how to attach uncertainties to a grid of data
# Uncertainties are the RMS in each channel, which get put on the data on each channel
apply_uncertainties = True
if apply_uncertainties:
  from uncertainties import unumpy

  print "Applying uncertainties to 252A (example of uncertainties coding)"

  if d["252A"].shape[1] != len(d["252A_rms"]):
    print "Mismatch in sizes: channels in uncertainties don't match data"
    exit(1)

  un_vals = np.zeros(d["252A"].shape)	# Prepare uncertainties array

  for i in range(un_vals.shape[1]):	# Fill uncertainties array
    un_vals[:, i] = d["252A_rms"][i]

  ant_252A_with_uncertainties = unumpy.uarray(d["252A"], un_vals)	# Make an array of data+uncertainties


