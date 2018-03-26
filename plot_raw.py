#!/usr/bin/env python
"""
# plot_raw.py

Plot the raw data by channel for 1 time step 
(the first in the file), unless specified
"""
import matplotlib as mpl
import os
import seaborn as sns
import tables as tb
from leda_cal.leda_cal import *
from leda_cal.skymodel import *
from leda_cal.dpflgr import *

sns.set_style('ticks')
sns.set_context("paper",font_scale=1.5)

def mid_range(h5, t):		# A range that is 3 standard deviations around the mean
  mn = 0
  rms = 0
  cols =[ h5.root.data.cols.ant252_x[t], h5.root.data.cols.ant252_x[t], h5.root.data.cols.ant252_x[t],
	h5.root.data.cols.ant252_x[t], h5.root.data.cols.ant252_x[t], h5.root.data.cols.ant252_x[t] ]
  for data in cols:
    mn += np.mean(data)
    rmsd = np.std(data)
    if rmsd > rms: rms = rmsd

  mn /= len(cols)
  return mn-3*rms, mn+3*rms


def quicklook(filename, save, noshow, t):
    h5 = tb.open_file(filename)

    if t >= len(h5.root.data.cols.ant252_x):
      print "Time step doesn't exist"
      exit(1)

    (bottom, top) = mid_range(h5, t)
    
    fig, ax = plt.subplots(figsize=(8, 6))

    plt.subplot(2,1,1)
    plt.plot(h5.root.data.cols.ant252_x[t], label="252A")
    plt.plot(h5.root.data.cols.ant254_x[t], label="254A")
    plt.plot(h5.root.data.cols.ant255_x[t], label="255A")
    plt.ylim(bottom, top)
    plt.legend(frameon=False)
    plt.ylabel("Data")
    plt.minorticks_on()
     
    plt.subplot(2,1,2)
    plt.plot(0, 0)
    plt.plot(h5.root.data.cols.ant252_y[t], label="252A")
    plt.plot(h5.root.data.cols.ant254_y[t], label="254B")
    plt.plot(h5.root.data.cols.ant255_y[t], label="255B")
    plt.ylim(bottom, top)
    plt.legend(frameon=False)
    
    plt.minorticks_on()
    
    
    plt.xlabel("Channel")
     
    plt.legend(frameon=False)
    plt.tight_layout()
    
    if save:
      plt.savefig("raw_"+os.path.basename(filename)[:-3]+".png")
    if not noshow:
      plt.show()

if __name__ == "__main__":
    import optparse, sys

    usage = '%prog [opts] filename_of_hdf5_observation'
    o = optparse.OptionParser()
    o.set_usage(usage)
    o.set_description(__doc__)
    o.add_option('--time', dest='time', default=0, type="int",
      help='The time step in the file. Default: 0')
    o.add_option('--save', dest='save', action='store_true', default=False,
      help="Save the plot to an image file, with filename the same as the h5 but png extension. Default: False.")
    o.add_option('--noshow', dest='noshow', action='store_true', default=False,
      help="Don't display the plot on screen. Useful for batch runs. Default: False.")

    opts, args = o.parse_args(sys.argv[1:])


    if len(args) != 1:
      o.print_help()
      exit(1)
    else: filename = args[0]

    quicklook(filename, opts.save, opts.noshow, opts.time)
 
