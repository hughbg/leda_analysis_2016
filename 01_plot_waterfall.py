#!/usr/bin/env python
"""
# 01_plot_waterfall.py

Plot data as calibrated waterfall plot.

"""
import os
import hickle
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import seaborn as sns
import tables as tb
import scipy.signal
from leda_cal.skymodel import *
from leda_cal.leda_cal import *
from leda_cal.dpflgr import *

sns.set_style('white')
sns.set_context("poster",font_scale=.75)

ovro_location = ('37.2397808', '-118.2816819', 1183.4839)
ovro = ephem.Observer(); (ovro.lat, ovro.lon, ovro.elev) = ovro_location
sun = ephem.Sun()

gal_center = ephem.FixedBody()  
gal_center._ra  =  '17 45 40.04'
gal_center._dec = '-29 00 28.1'
gal_center.name = "Galactic Center"


def quicklook(filename, save, dump, flag, no_show, all_lsts):
    h5 = tb.open_file(filename)

    T_ant = apply_calibration(h5)
    f_leda = T_ant['f']
    
    ant_ids = ['252', '254', '255']
      
    print("Plotting...")
    fig, axes = plt.subplots(figsize=(12, 12))
    #plt.suptitle(h5.filename)
    
    lst_stamps = T_ant['lst']
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
       padding = np.delete(padding, unusable_lsts, axis=0)
       lst_stamps = np.delete(lst_stamps, unusable_lsts, axis=0)
       utc_stamps = np.delete(utc_stamps, unusable_lsts, axis=0)
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

    for p in [ 0, 1 ]:

      for ii, key in enumerate(ant_ids):
          ax = fig.add_subplot(2, 3, 3*p+ii+1)
        
	  if p == 0: ant = key+"A"
          else: ant = key+"B"

          T_flagged = T_ant[ant]
          if not all_lsts:  T_flagged = np.delete(T_flagged, unusable_lsts, axis=0)
 
          if flag:
            T_flagged = rfi_flag(T_flagged, thr_f=0.2, thr_t=0.2, rho=1.5,
               bp_window_f=16, bp_window_t=16,
               max_frac_f=0.5, max_frac_t=0.5)


          
          if dump: dump_data[ant] = T_flagged
          if flag: 
            total = T_flagged.shape[0]*T_flagged.shape[1]
            num_in = np.ma.MaskedArray.count(T_flagged)
            print ant, ( "%.1f%%" % (100*(total-num_in)/total) ), "flagged.", "Count:", total-num_in
        
          # Add the stripe onto the right edge of the data and adjust the extent of the x-axis (frequency) to cover the stripe.
          if all_lsts: T_flagged_plot = np.ma.concatenate((T_flagged, padding), axis=1)
          else: T_flagged_plot = T_flagged

          axim = plt.subplot(2, 3, 3*p+ii+1)
          axim.set_yticks(yloc)
          axim.set_yticklabels(ylabel)
	  axim.tick_params(axis='y', pad=2)

          im = plt.imshow(T_flagged_plot, # / np.median(xx, axis=0), 
                   cmap='viridis', aspect='auto',
                   interpolation='nearest',
                   clim=(1000, 10000),
                   extent=(xlims[0], new_x_high, ylims[1], ylims[0])
                   )

          plt.title(ant)
 	  if p == 1: plt.xlabel("Frequency [MHz]")
          #ax.yaxis_date()
          #ax.yaxis.set_major_formatter(hfmt)
          #

    plt.subplot(2,3,1)
    plt.ylabel("LST [hr]")
    plt.subplot(2,3,4)
    plt.ylabel("LST [hr]")

    
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
      dump_data["options"] = "Flag="+str(flag)+" All LSTSs="+str(all_lsts)
      hickle.dump(dump_data, os.path.basename(filename)[:-3]+".hkl")

if __name__ == "__main__":
    import optparse, sys

    usage = '%prog [opts] filename_of_hdf5_observation'
    o = optparse.OptionParser()
    o.set_usage(usage)
    o.set_description(__doc__)
    o.add_option('--flag', dest='flag', action='store_true', default=False,
      help='Apply flagging. Default: False')
    o.add_option('--all_lsts', dest='all_lsts', action='store_true', default=False,
      help='Include all LSTs, not just when Galaxy and Sun are down. A day/night stripe is printed on the right. Default: False.')
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

    quicklook(filename, opts.save, opts.dump, opts.flag, opts.no_show, opts.all_lsts)
