#!/usr/bin/env python
"""
# 01_plot_waterfall.py

Plot data as calibrated waterfall plot.

"""
import os
import hickle
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import seaborn as sns
import tables as tb
import scipy.signal
from scipy.stats import scoreatpercentile as percentile
from leda_cal.skymodel import *
from leda_cal.leda_cal import *
from leda_cal.dpflgr import *
from leda_cal.useful import add_uncertainties
import leda_cal.params

sns.set_style('white')
sns.set_context("poster",font_scale=.75)

ovro_location = ('37.2397808', '-118.2816819', 1183.4839)
ovro = ephem.Observer(); (ovro.lat, ovro.lon, ovro.elev) = ovro_location
sun = ephem.Sun()

gal_center = ephem.FixedBody()  
gal_center._ra  =  '17 45 40.04'
gal_center._dec = '-29 00 28.1'
gal_center.name = "Galactic Center"


def quicklook(filename, save, dump, flag, merge, flatten, no_show, all_lsts):
    h5 = tb.open_file(filename)

    T_ant = apply_calibration(h5)
    f_leda = T_ant['f']
    
    ant_ids = ['252', '254', '255']
      
    print("Plotting...")
    fig = plt.figure(figsize=(12, 12))
    #plt.suptitle(h5.filename)
    
    lst_stamps = T_ant['lst']
    if len(lst_stamps) == 0:
      print "No LSTS in file"
      exit(1)

    # Report discontinuities in time
    for i in range(1,len(lst_stamps)):
      if lst_stamps[i]-lst_stamps[i-1] > 1/60.0:	# 1 minute
        print "Discontinuity at LST", lst_stamps[i], (lst_stamps[i]-lst_stamps[i-1])*60*60, "seconds"

    utc_stamps = T_ant['utc']
    xlims = (f_leda[0], f_leda[-1])
    #ylims = mdates.date2num((T_ant['utc'][0], T_ant['utc'][-1]))
    #hfmt = mdates.DateFormatter('%m/%d %H:%M')
    ylims = (T_ant['lst'][0], T_ant['lst'][-1])

    # Work out altitude of Gal center and Sun. Use whichever is highest
    # and put that in the padding, which is the stripe.
    unusable_lsts = []
    pad_length = 70
    padding = np.zeros((len(lst_stamps), pad_length))
    for i, d in enumerate(utc_stamps):
      ovro.date = d
      sun.compute(ovro)
      gal_center.compute(ovro)
      if sun.alt > -15*np.pi/180 or gal_center.alt > -15*np.pi/180:
        padding[i, :] = 10000
        unusable_lsts.append(i)
      else: 
        padding[i, :] = 1000
 
    # Delete sun up LSTS
    if not all_lsts:
       print "Cutting out times when sun/galaxy up"
       padding = np.delete(padding, unusable_lsts, axis=0)
       lst_stamps = np.delete(lst_stamps, unusable_lsts, axis=0)
       utc_stamps = np.delete(utc_stamps, unusable_lsts, axis=0)
       if len(lst_stamps) == 0:
         print "No LSTs available at night time (use --all_lsts to see all)"
         exit(1)
       ylims = ( lst_stamps[0], lst_stamps[-1] )
       print len(lst_stamps), "usable LSTs"
    else: print "Using all LSTs"
    if len(lst_stamps) == 0:
      print "There is no data to display (number of LSTs is 0)"
      exit(1) 

    yloc = []
    ylabel = []
    for i in range(0, len(lst_stamps), len(lst_stamps)/7):
      yloc.append(lst_stamps[i]), ylabel.append(("%.1f" % lst_stamps[i]))
    if all_lsts: new_x_high = xlims[1]+pad_length*(xlims[1]-xlims[0])/len(f_leda)
    else: new_x_high = xlims[1]

    dump_data = {}
    
    if flag and merge:
        # If we are going to merge the flags across antennas, we need to flag them all now
        for p in (0, 1):
            for ii,key in enumerate(ant_ids):
                ant = key+("B" if p else "A")
                T_flagged = T_ant[ant]
                if not all_lsts:
                    T_flagged = np.delete(T_flagged, unusable_lsts, axis=0)
                new_mask = rfi_flag(T_flagged, freqs=f_leda).mask
                try:
                    merged_mask |= new_mask
                except NameError:
                    merged_mask = new_mask
        
    for p in [ 0, 1 ]:

      for ii, key in enumerate(ant_ids):
          if p == 0 and ii == 0:
              ax = fig.add_subplot(2, 3, 3*p+ii+1)
              origAX = ax
          else:
              ax = fig.add_subplot(2, 3, 3*p+ii+1, sharex=origAX, sharey=origAX)
              
	  if p == 0: ant = key+"A"
          else: ant = key+"B"

          T_flagged = T_ant[ant]
          if not all_lsts:  T_flagged = np.delete(T_flagged, unusable_lsts, axis=0)

          print "Max", np.max(T_flagged), "Min", np.min(T_flagged)

          if flag:
            if merge:
                ## Already done
                T_flagged = np.ma.array(T_flagged, mask=merged_mask)
            else:
                ## Need to do it now - there's probably a way to deal with 
                ## this all in one pass
                T_flagged = rfi_flag(T_flagged, freqs=f_leda)
            print "After flagging", "Max", np.ma.max(T_flagged), "Min", np.ma.min(T_flagged)

          if dump: 
            dump_data[ant] = T_flagged
            dump_data[ant+"_rms"] = add_uncertainties(T_flagged)
	    av = np.ma.average(T_flagged,axis=0)
	    weighted = av/dump_data[ant+"_rms"]**2
            dump_data[ant+"_weighted"]=weighted       
            
          
          if flag: 
            total = T_flagged.shape[0]*T_flagged.shape[1]
            num_in = np.ma.MaskedArray.count(T_flagged)
            print ant, ( "%.1f%%" % (100*(total-num_in)/total) ), "flagged.", "Count:", total-num_in
        
          # Add the stripe onto the right edge of the data and adjust the extent of the x-axis (frequency) to cover the stripe.
          if all_lsts: T_flagged_plot = np.ma.concatenate((T_flagged, padding), axis=1)
          else: T_flagged_plot = T_flagged

          #if p == 0 and ii == 0:
          #    axim = origAX
          #else:
          #    axim = fig.add_subplot(2, 3, 3*p+ii+1, sharex=origAX, sharey=origAX)
          ax.set_yticks(yloc)
          ax.set_yticklabels(ylabel)
          ax.tick_params(axis='y', pad=2)
          
          if flatten:
             abp = np.ma.median(T_flagged_plot, axis=0)
             abp /= np.ma.median(abp)
             T_flagged_plot /= abp
             try:
                 clim = (percentile(T_flagged_plot.compressed(), 5), percentile(T_flagged_plot.compressed(), 95))
             except AttributeError:
                 clim = (percentile(T_flagged_plot, 5), percentile(T_flagged_plot, 95))

          else:
             clim = (1000, 10000)

          im = ax.imshow(T_flagged_plot, # / np.median(xx, axis=0), 
                   cmap='jet', aspect='auto',
                   interpolation='nearest',
                   clim=clim,
                   extent=(xlims[0], new_x_high, ylims[1], ylims[0])
                   )

          ax.set_title(ant)
          if p == 1: ax.set_xlabel("Frequency [MHz]")
          if ii == 0: ax.set_ylabel("LST [hr]")
          #ax.yaxis_date()
          #ax.yaxis.set_major_formatter(hfmt)
          #

    if not flatten:
        fig.subplots_adjust(left=0.07)
        fig.subplots_adjust(right=0.875)
        cbar_ax = fig.add_axes([0.9, 0.125, 0.025, 0.75])
        cbar = fig.colorbar(im, cax=cbar_ax)
        
        #plt.subplot(2,3,3)
        #cbar = plt.colorbar()
        cbar.set_label("Temperature [K]")
        cbar.ax.tick_params(axis='y', pad=2) 
        #plt.tight_layout()
    
    if save:
      plt.savefig(os.path.basename(filename)[:-3]+".png")
    if not no_show:
      plt.show()
    
    if dump:
      dump_data["lsts"] = lst_stamps
      dump_data["utcs"] = np.array([str(pytime) for pytime in utc_stamps])
      dump_data["frequencies"] = f_leda
      dump_data["options"] = "Flag="+str(flag)+" Merge="+str(merge)+" Flatten="+str(flatten)+" All LSTSs="+str(all_lsts)
      hickle.dump(dump_data, os.path.basename(filename)[:-3]+".hkl")

if __name__ == "__main__":
    import optparse, sys

    usage = '%prog [opts] filename_of_hdf5_observation'
    o = optparse.OptionParser()
    o.set_usage(usage)
    o.set_description(__doc__)
    o.add_option('--flatten', dest='flatten', action='store_true', default=False,
      help='Apply a crude bandpass derived from the data. Default: False')
    o.add_option('--flag', dest='flag', action='store_true', default=False,
      help='Apply flagging. Default: False')
    o.add_option('--merge', dest='merge', action='store_true', default=False,
      help='Merge all flags. Default: False')
    o.add_option('--all_lsts', dest='all_lsts', action='store_true', default=False,
      help='Include all LSTs, not just when Galaxy and Sun are down. A day/night stripe is printed on the right. Default: False.')
    o.add_option('--median', dest='median', action='store_true', default=False,
      help='Use a median filter in the sum threshold flagging. Overrides setting in params. Default: False.')
    o.add_option('--dump', dest='dump', action='store_true', default=False,
      help='Dump the data to a hickle file, with filename the same as the h5 but hkl extension. Default: False.')
    o.add_option('--save', dest='save', action='store_true', default=False,
      help="Save the plot to an image file, with filename the same as the h5 but png extension. Default: False.")
    o.add_option('--no_show', dest='no_show', action='store_true', default=False,
      help="Don't display the plot on screen. Useful for batch runs. Default: False.")

    opts, args = o.parse_args(sys.argv[1:])
    
    if len(args) != 1:
      o.print_help()
      exit(1)
    else: filename = args[0]

    params.median = opts.median
    quicklook(filename, opts.save, opts.dump, opts.flag, opts.merge, opts.flatten, opts.no_show, opts.all_lsts)
