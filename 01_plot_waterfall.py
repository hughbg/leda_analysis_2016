#!/usr/bin/env python
"""
# 01_plot_waterfall.py

Plot data as calibrated waterfall plot.

"""
import os
import hickle
import numpy as np
import sys
if "--no_show" in sys.argv: 
  from matplotlib import use as muse; muse('Agg')
from matplotlib import cm
import seaborn as sns
import tables as tb
import scipy.signal
from scipy.stats import scoreatpercentile as percentile
from leda_cal.skymodel import *
from leda_cal.leda_cal import *
from leda_cal.dpflgr import *
from leda_cal.useful import add_uncertainties, ensure_mask
from leda_cal import robust
from leda_cal.git import get_repo_fingerprint
from leda_cal.params import params
from leda_cal import lst_timing

sns.set_style('white')
sns.set_context("poster",font_scale=.75)

REQUIRED = 2400		# channels. Pad if necessary.

def biggest_gap(times):	# times must be sorted
    gap = -1
    where_gap = (-1, -1)
    for i in range(1, len(times)):
      this_gap = times[i]-times[i-1]
      if this_gap > gap:
        gap = this_gap
        where_gap = (times[i], times[i-1])
    return where_gap

def pad_data(data):


  if data.shape[1] == REQUIRED: return data

  missing = REQUIRED-data.shape[1]
  if missing > 0:
    pad = np.ma.masked_all((data.shape[0], missing))
    new_data = np.ma.concatenate((ensure_mask(data), pad), axis=1)

  return new_data

def pad_frequencies(frequencies):

  last_index = len(frequencies)-1
  last_val = frequencies[-1]
  chan_width = (frequencies[-1]-frequencies[0])/len(frequencies-1)
  frequencies = np.append(frequencies, np.zeros(REQUIRED-len(frequencies)))
  for i in range(last_index+1, len(frequencies)):
     frequencies[i] = last_val+(i-last_index)*chan_width

  return frequencies


def quicklook(filename, save, dump, flag, merge, flatten, no_show, all_lsts, new_cal, sky=False, lfsm=False, emp=False):
    h5 = tb.open_file(filename)
 
    if new_cal: T_ant = apply_new_calibration(h5)
    else: T_ant = apply_calibration(h5)
    f_leda = T_ant['f']
    
    ant_ids = ['252', '254', '255']
    
    print("Plotting...")
    fig = plt.figure(figsize=(20,20))
    #plt.suptitle(h5.filename)
    
    lst_stamps = T_ant['lst']
    indexes = np.arange(len(lst_stamps), dtype=np.int)

    if len(lst_stamps) == 0:
        raise RuntimeError("No LSTs in file")

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
    pad_length = 70
    padding = np.full((len(lst_stamps), pad_length), 10000)
    timing = lst_timing.LST_Timing(lst_stamps, utc_stamps)	
    border_bottom, night_bottom, night_top, border_top = timing.calc_night()
    padding[night_bottom:night_top, :] = 1000

    #for ant in ant_ids:
    #  lst_stamps, T_ant[ant+"A"] = timing.align(T_ant[ant+"A"])
    #  lst_stamps, T_ant[ant+"B"] = timing.align(T_ant[ant+"B"])

    if night_bottom: print "Night", lst_stamps[night_bottom], "-", lst_stamps[night_top-1]
    else: print "Night 0 - 0"

    # Use night only
    if not all_lsts:
      if not border_top:
        raise RuntimeError("No LSTs available at night time (use --all_lsts to see all)")
      lst_stamps = lst_stamps[night_bottom:night_top]
      utc_stamps = utc_stamps[night_bottom:night_top]
      indexes = indexes[night_bottom:night_top]
      padding = padding[night_bottom:night_top]
      ylims = ( lst_stamps[0], lst_stamps[-1] )
      print len(lst_stamps), "usable LSTs"
    else:
      print "Using all LSTs"
   

    if len(lst_stamps) == 0:
        raise RuntimeError("There are no data to display (number of LSTs is 0)")
        
    yloc = []
    ylabel = []
    try:
      for i in range(0, len(lst_stamps), len(lst_stamps)/7):
        yloc.append(lst_stamps[i]), ylabel.append(("%.1f" % lst_stamps[i]))
    except: 
      yloc.append(lst_stamps[0]), ylabel.append(("%.1f" % lst_stamps[0]))
      yloc.append(lst_stamps[-1]), ylabel.append(("%.1f" % lst_stamps[-1]))
    if all_lsts:
        new_x_high = xlims[1]+pad_length*(xlims[1]-xlims[0])/len(f_leda)
    else:
        new_x_high = xlims[1]
        
    dump_data = {}
    
    if sky:
        if lfsm and emp:
            smdl = SkyModelLFSMEmp
            smlbl = 'LFSM+Emp'
        elif lfsm and not emp:
            smdl = SkyModelLFSM
            smlbl = 'LFSM'
        elif not lfsm and emp:
            smdl = SkyModelGSMEmp
            smlbl = 'GSM+Emp'
        else:
            smdl = SkyModelGSM        
            smlbl = 'GSM'
        sy = smdl(pol='y')
        sx = smdl(pol='x')
        T_y_asm = sy.generate_tsky(lst_stamps, f_leda*1e6)
        T_x_asm = sx.generate_tsky(lst_stamps, f_leda*1e6)

    dtv_times = {}
    if flag and merge:
        # If we are going to merge the flags across antennas, we need to flag them all now
        for p in (0, 1):
            for ii,key in enumerate(ant_ids):
                ant = key+("B" if p else "A")
                T_flagged = T_ant[ant]
                if not all_lsts: 
		  # Do flagging with a border around the data in time
		  flagged, dtv_tms = rfi_flag(T_flagged[border_bottom:border_top], freqs=f_leda)
		  new_mask = flagged.mask	
		  
                  new_mask = new_mask[night_bottom-border_bottom:night_top-border_bottom]		# remove border
   	          dtv_times[ant] = [ x-(night_bottom-border_bottom) for x in dtv_tms ]
		  dtv_times[ant] = [ x for x in dtv_times[ant] if 0 <= x and x < night_top-night_bottom ] 		# Some might have been in the border
		else:
		  flagged, tms = rfi_flag(T_flagged, freqs=f_leda)
                  new_mask = flagged.mask
		  dtv_times[ant] = tms
		  print ant, "Biggest DTV gap", lst_stamps[biggest_gap(dtv_tms)[1]], "-", lst_stamps[biggest_gap(dtv_tms)[0]], "waterfall"
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
                
            if p == 0:
                ant = key+"A"
            else:
                ant = key+"B"
                
            T_flagged = T_ant[ant]
            if not all_lsts:
                T_flagged = T_flagged[night_bottom:night_top]
                

            print "Max", np.max(T_flagged), "Min", np.min(T_flagged)
            
	    dtv_times[ant] = []
            if flag:
                if merge:
                    ## Already done
                    T_flagged = np.ma.array(T_flagged, mask=merged_mask)
                else:
                    ## Need to do it now - there's probably a way to deal with 
                    ## this all in one pass
		    if not all_lsts:
		      flagged, dtv_tms = rfi_flag(T_ant[ant][border_bottom:border_top], freqs=f_leda)
		      T_flagged = flagged
		      T_flagged = T_flagged[night_bottom-border_bottom:night_top-border_bottom]	# Remove border
		      dtv_times[ant] = [ x-(night_bottom-border_bottom) for x in dtv_tms ]
		      dtv_times[ant] = [ x for x in dtv_times[ant] if 0 <= x and x < night_top-night_bottom ] 		# Some might have been in the border
		    else: 
                      T_flagged, dtv_tms = rfi_flag(T_flagged, freqs=f_leda)
		      dtv_times[ant] = dtv_tms
		      print ant, "Biggest DTV gap", lst_stamps[biggest_gap(dtv_tms)[1]], "-", lst_stamps[biggest_gap(dtv_tms)[0]], "waterfall"
                print "After flagging", "Max", np.ma.max(T_flagged), "Min", np.ma.min(T_flagged)
                
            try:
                T_asm = T_y_asm if p == 0 else T_x_asm
                scale_offset_asm = robust.mean(T_asm / T_flagged)
                T_flagged = T_flagged - T_asm / scale_offset_asm
            except NameError:
                pass

	    T_flagged = pad_data(T_flagged)		# Up to 2400 channels

            if dump: 
                dump_data[ant] = T_flagged
                dump_data[ant+"_rms"] = add_uncertainties(T_flagged)
                av = np.ma.average(T_flagged,axis=0)
                weighted = av/dump_data[ant+"_rms"]**2
                dump_data[ant+"_weighted"] = weighted    
	        dump_data[ant+"_dtv_times"] = np.array(dtv_times[ant], dtype=np.int) 
                
            if flag:
                total = T_flagged.shape[0]*T_flagged.shape[1]
                num_in = np.ma.MaskedArray.count(T_flagged)
                print ant, ( "%.1f%%" % (100*float(total-num_in)/total) ), "flagged.", "Count:", total-num_in
                
            # Add the stripe onto the right edge of the data and adjust the extent of the x-axis (frequency) to cover the stripe.
            if all_lsts:
                T_flagged_plot = np.ma.concatenate((T_flagged, padding), axis=1)
            else:
                T_flagged_plot = T_flagged
                
            ax.set_yticks(yloc)
            ax.set_yticklabels(ylabel)
            ax.tick_params(axis='y', pad=2)
            
            if flatten:
                if type(T_flagged_plot) is np.ma.core.MaskedArray:
                    abp = np.ma.median(T_flagged_plot.data, axis=0)
                else:
                    abp = np.ma.median(T_flagged_plot, axis=0)
                abp /= np.ma.median(abp)
                T_flagged_plot /= abp
                try:
                    clim = (percentile(T_flagged_plot.compressed(), 5), percentile(T_flagged_plot.compressed(), 95))
                except AttributeError:
                    clim = (percentile(T_flagged_plot, 5), percentile(T_flagged_plot, 95))

            elif sky:
                clim = (-250, 500)
            else:
                clim = (1000, 10000)

            im = ax.imshow(T_flagged_plot, # / np.median(xx, axis=0), 
                           cmap="viridis", aspect='auto',
                           interpolation='nearest',
                           clim=clim,
                           extent=(xlims[0], new_x_high, ylims[1], ylims[0]))
            
            ax.set_title(ant)
            if p == 1:
                ax.set_xlabel("Frequency [MHz]")
            if ii == 0:
                ax.set_ylabel("LST [hr]")
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
        if sky:
            cbar.set_label("Temperature - %s [K]" % smlbl)
        else:
            cbar.set_label("Temperature [K]")
        cbar.ax.tick_params(axis='y', pad=2) 
        #plt.tight_layout()
    
    plt.text(0.005, 0.005, get_repo_fingerprint(), transform=fig.transFigure, size=8)

    if save:
        plt.savefig(os.path.basename(filename)[:-3]+".png")
    if not no_show:
        plt.show()
        
    if dump:
        dump_data["lsts"] = lst_stamps
        dump_data["utcs"] = np.array([str(pytime) for pytime in utc_stamps])
        dump_data["indexes"] = indexes
        dump_data["frequencies"] = pad_frequencies(f_leda)
        dump_data["options"] = "Flag="+str(flag) \
			       + " Filename="+filename \
			       + " New cal="+str(new_cal) \
                               + " Merge="+str(merge) \
                               + " Flatten="+str(flatten) \
                               + " All LSTs="+str(all_lsts) \
                               + " Sky Model Substract="+str(sky) \
                               + " Use LFSM="+str(lfsm) \
                               + " Apply empirical gain correction="+str(emp)
        dump_data["fingerprint"] = get_repo_fingerprint()
        import json
        def jdefault(o): return o.__dict__
        dump_data["params"] = json.dumps(params, default=jdefault)

        hickle.dump(dump_data, os.path.basename(filename)[:-3]+".hkl")


if __name__ == "__main__":
    import optparse

    usage = '%prog [opts] filename_of_hdf5_observation'
    o = optparse.OptionParser()
    o.set_usage(usage)
    o.set_description(__doc__)
    o.add_option('--sky', dest='sky', action='store_true', default=False,
                 help='Subtract off a sky model. Default: False')
    o.add_option('--lfsm', dest='lfsm', action='store_true', default=False,
                 help='Use the LFSM instead of the GSM')
    o.add_option('--empirical', dest='emp', action='store_true', default=False,
                 help='Apply an empirical corretion to the dipole gain pattern model')
    o.add_option('--flatten', dest='flatten', action='store_true', default=False,
                 help='Apply a crude bandpass derived from the data. Default: False')
    o.add_option('--flag', dest='flag', action='store_true', default=False,
                 help='Apply flagging. Default: False')
    o.add_option('--merge', dest='merge', action='store_true', default=False,
                 help='Merge all flags. Default: False')
    o.add_option('--new_cal', dest='new_cal', action='store_true', default=False,
                 help='Use the 2018 calibration files and method. Default: False.')   
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
    else:
        filename = args[0]
        
    params.median = opts.median
    quicklook(filename, opts.save, opts.dump, opts.flag, opts.merge, opts.flatten, opts.no_show, opts.all_lsts, opts.new_cal,
              sky=opts.sky, lfsm=opts.lfsm, emp=opts.emp)
