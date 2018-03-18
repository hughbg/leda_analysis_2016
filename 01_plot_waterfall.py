#!/usr/bin/env python
"""
# 01_plot_waterfall.py

Plot data as calibrated waterfall plot.

"""
import os
import seaborn as sns
import tables as tb
from leda_cal.skymodel import *
from leda_cal.leda_cal import *

sns.set_style('white')
sns.set_context("poster",font_scale=.75)

ovro_location = ('37.2397808', '-118.2816819', 1183.4839)
ovro = ephem.Observer(); (ovro.lat, ovro.lon, ovro.elev) = ovro_location
sun = ephem.Sun()

def quicklook(filename, pol, save, dump):
    h5 = tb.open_file(filename)

    T_ant = apply_calibration(h5)
    f_leda = T_ant['f']
    
    ant_ids = ['252'+pol, '254'+pol, '255'+pol]
    pol_id  = 'y'
      
    print("Plotting...")
    fig, axes = plt.subplots(figsize=(12, 6), nrows=2, ncols=2)
    #plt.suptitle(h5.filename)
    
    lst_stamps = T_ant['lst']
    utc_stamps = T_ant['utc']
    xlims = (f_leda[0], f_leda[-1])
    #ylims = mdates.date2num((T_ant['utc'][0], T_ant['utc'][-1]))
    #hfmt = mdates.DateFormatter('%m/%d %H:%M')
    ylims = (T_ant['lst'][0], T_ant['lst'][-1])
    for ii, key in enumerate(ant_ids):
        ax = fig.add_subplot(1, 3, ii+1)
        
        T_flagged = T_ant[key]
        #T_flagged = rfi_flag(T_ant[key], thr_f=0.07, thr_t=0.07, rho=1.5,
        #         bp_window_f=8, bp_window_t=8, 
        #         max_frac_f=0.3, max_frac_t=0.3)
        
        T_ant[key] = T_flagged
        
	if dump: np.save(os.path.basename(filename)[:-3]+"_"+key, T_flagged)
        max_alt = 0
        pad_length = 70
        padding = np.zeros((T_flagged.shape[0], pad_length))
        for i, d in enumerate(utc_stamps):
          ovro.date = d
          sun.compute(ovro)
          padding[i, :] = sun.alt
          if sun.alt > 0 and sun.alt > max_alt: max_alt = abs(sun.alt)

        for i in range(len(utc_stamps)):
          padding[i, :] = 1000+((padding[i, 0]+max_alt)/(2*max_alt))*(10000-1000)
          #print d, lst_stamps[i], ovro.sidereal_time(), float(sun.alt)*180/np.pi, 1000+((sun.alt+max_alt)/(2*max_alt))*(10000-1000)
          
        T_flagged = np.concatenate((T_flagged, padding), axis=1)
        new_x_high = xlims[1]+pad_length*(xlims[1]-xlims[0])/T_flagged.shape[1]
        im = plt.imshow(T_flagged, # / np.median(xx, axis=0), 
                   cmap='viridis', aspect='auto',
                   interpolation='nearest',
                   clim=(1000, 10000),
                   extent=(xlims[0], new_x_high, ylims[1], ylims[0])
                   )
        plt.title(ant_ids[ii])
 	plt.xlabel("Frequency [MHz]")
        #ax.yaxis_date()
        #ax.yaxis.set_major_formatter(hfmt)
        #

    
    plt.subplot(1,3,1)
    plt.ylabel("LST [hr]")
    
    fig.subplots_adjust(left=0.06)
    fig.subplots_adjust(right=0.875)
    cbar_ax = fig.add_axes([0.9, 0.125, 0.025, 0.75])
    cbar = fig.colorbar(im, cax=cbar_ax)
    
    plt.subplot(1,3,3)
    #cbar = plt.colorbar()
    cbar.set_label("Temperature [K]")
    
    #plt.tight_layout()
    
    if save:
      plt.savefig(os.path.basename(filename)[:-3]+"_"+pol+".png")
    else: 
      plt.show()
    
    

if __name__ == "__main__":
    import optparse, sys

    usage = '%prog [opts] filename_of_hdf5_observation'
    o = optparse.OptionParser()
    o.set_usage(usage)
    o.set_description(__doc__)
    o.add_option('--pol', dest='pol', default='A',
      help='Polarization A or B. Default: A')
    o.add_option('--dump', dest='dump', action='store_true', default=False,
      help='Dump the data to a file, with filename the same as the h5 but npy extension. Default: False')
    o.add_option('--save', dest='save', action='store_true', default=False,
      help='Save the plot to an image file, with filename the same as the h5 but png extension. Default: False')
    opts, args = o.parse_args(sys.argv[1:])
    
    if len(args) != 1:
      o.print_help()
      exit(1)
    else: filename = args[0]
    if opts.pol not in [ "A", "B" ]:
      print "Invalid pol"
      o.print_help()
      exit(1)

    quicklook(filename, opts.pol, opts.save, opts.dump)
