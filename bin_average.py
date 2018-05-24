#!/usr/bin/env python

#  Usage like: python bin_average.py outriggers_2018-04-02_07H52M55S.hkl 252 9 30
# Creates a hickle file with averages, in the dav directory. Read them with dav/hkl_report.py

import hickle, sys, numpy as np
import matplotlib.pyplot as plt
import math
import bottleneck as bn
import scipy.optimize


if len(sys.argv) != 5:
  print "Usage:", sys.argv[0], "hkl_file  antenna  LST_start av_length"
  print "	hkl_file:	dumped from 01_plot_waterfall.py"
  print "	antenna:	like 252, without the A/B"
  print "	LST_start:	an LST that exists in the data"
  print "	av_length:	time in hours to bin by, e.g. 1"
  exit(1)

def scrunch(a): 		# Do this is in a better numpy way. useful.rebin leaves off end values 
  num_in_bin = 42		# 1MHz


  if a.ndim == 1:
    start = a.shape[0]-(a.shape[0]/num_in_bin)*num_in_bin
    result = np.zeros(a.shape[0]/num_in_bin)
    for i in range(start, a.shape[0], num_in_bin):  
      result[i/num_in_bin] = np.ma.average(a[i: i+num_in_bin])

  else:
    start = a.shape[1]-(a.shape[1]/num_in_bin)*num_in_bin
    result = np.zeros((a.shape[0], a.shape[1]/num_in_bin))
    for i in range(start, a.shape[1], num_in_bin):
      result[:, i/num_in_bin] = np.ma.average(a[:, i: i+num_in_bin], axis=1)

  return result

# Find the spectral index just by curve fitting
def spectral_index(freq, temp):
  def S(x, C, si):
    return C*x**si

  # Don't want masked values or NaNs or zeros
  nan_temp = np.ma.filled(temp, np.nan)
  f = freq[np.logical_or(np.isnan(nan_temp), (nan_temp!=0))]
  s = temp[np.logical_or(np.isnan(nan_temp), (nan_temp!=0))]

  try:
    popt, pcov = scipy.optimize.curve_fit(S, f, s)
  except: popt = ( 0, 0 )

  return popt[1]

def load_hkl(filename):
  return hickle.load(filename)
d = load_hkl(sys.argv[1])


print "LST", d["lsts"][0], "to", d["lsts"][-1]

dipoles = ['A','B']
min_lst = float(sys.argv[3]) # in hours
dT = float(sys.argv[4])# in hours
vecTbins = np.arange(min_lst,max(d["lsts"]),dT)
use_bins = []
for ind in range(len(vecTbins)-1):
  lst_ind = [i for i in range(len(d["lsts"])) if d["lsts"][i]>=vecTbins[ind] and d["lsts"][i]<vecTbins[ind+1]]
 
  if len(lst_ind) > 0: bin_ok = True
  else: bin_ok = False
  for i in lst_ind[1:]:
    if d["lsts"][i]-d["lsts"][i-1] > 1/60.0:	# 1 minute
      print "Discontinuity at LST", d["lsts"][i], (d["lsts"][i]-d["lsts"][i-1])*60*60, "seconds"
      bin_ok = False
  use_bins.append(bin_ok)

Time_bins = []
bp_window_t = 8
av_data_dict = {}
av_data_dict["Bins"] = [ {} for i in range(len(vecTbins)-1) ]
av_data_dict["Freqs"] = scrunch(d["frequencies"])
fig = plt.figure() 
	


flagsp = 1
for indD in range(len(dipoles)):

	ant = sys.argv[2]+dipoles[indD]		#"255A"

#print "file", sys.argv[1][0:len(sys.argv[1])-4]
#tlst = d["lsts"]# time in hours
	
#print "Length Freq",len(d["frequencies"])
#print "D", d[ant][1000][100]
#print "min T", min(d["lsts"])
#print "max T", max(d["lsts"])


#print "T_bins", vecTbins

	freq_scrunched = scrunch(d["frequencies"])

	data_averaged =  np.ma.zeros(len(freq_scrunched))
	rms_av = np.zeros(len(freq_scrunched))
	
	for ind in range(len(vecTbins)-1):

 		if not use_bins[ind]: continue

		
		indexT = [i for i in range(len(d["lsts"])) if d["lsts"][i]>=vecTbins[ind] and d["lsts"][i]<vecTbins[ind+1]]

		av_data_dict["Bins"][ind]["Time"] = (d["lsts"][indexT[0]], d["lsts"][indexT[-1]], len(indexT))

		d_scrunched = scrunch(d[ant][indexT])					# Averages over 1MHz freq bins.  Masked entries are ignored
		data_averaged = np.ma.average(d_scrunched, axis=0)  			# averaging over time
											# seems like if channel is all masked it becomes 0
  		num_in = np.ma.MaskedArray.count(data_averaged)
  		if num_in != len(data_averaged):
    		  print "ERROR: There are masked values in averaged array"
	          exit(1)
		if len(data_averaged[np.isnan(data_averaged)]) != 0:
   		  print "ERROR: There are NaN values in averaged array"
		  exit(1)

 		av_data_dict["Bins"][ind][ant] = data_averaged

		flat = bn.move_nanmean(d_scrunched, bp_window_t, axis=0)
		flat = np.roll(flat, -bp_window_t/2+1, axis=0)
		for indnu in range(d_scrunched.shape[1]):
			#print "indnu", d[ant][indexT,indnu].shape
			flat_channel = d_scrunched[:, indnu]-flat[ :, indnu]
			rms_av[indnu] = np.std(flat_channel[bp_window_t:])
			
	#		print "rms",rms_av[ind][indnu]
	#	plt.subplot(len(vecTbins)-1,2,flagsp-2*(len(vecTbins)-1)*indD+indD)
	#	plt.subplot(1,2,flagsp-2*(len(vecTbins)-1)*indD+indD)
		
		av_data_dict["Bins"][ind][ant+"_RMS"] = rms_av
		av_data_dict["Bins"][ind][ant+"_SI"] = spectral_index(freq_scrunched[19:], data_averaged[19:])
		
		#plt.subplot(2,2,flagsp)
	
		#flagsp=flagsp+1
		#flagsp=flagsp+2
		#plt.plot(d["frequencies"],data_averaged[ind][:])
		#plt.errorbar(d["frequencies"],data_averaged[ind][:], yerr = rms_av[ind][:])	
#		if ind==0: 
		#plt.ylabel('T [K]/1000 K')
		#plt.ylim(0,8000)
		#plt.yticks([0,2500,5000,7500],[0,2.5,5,7.5])
		#plt.xlim(30,90)
		#plt.xlabel('Freq. [MHz]')
		#plt.title("Dipole: "+ dipoles[indD]+ "   T = " +str(round(Time_bins[ind],2))+"[h]")
		
		#axim = plt.subplot(2, 2, flagsp+2)
 		#flagsp=flagsp+1
          	#plt.xlabel('Freq. [MHz]')
		#plt.ylabel('Time [h, lst]')
		#axim.set_yticks(yloc)
          	#axim.set_yticklabels(ylabel)
	  	#axim.tick_params(axis='y', pad=2)

          	#im = plt.imshow(d[ant][indexT][:], # / np.median(xx, axis=0), 
                #   cmap='viridis', aspect='auto',
                #   interpolation='nearest',
                #   clim=(1000, 7500),
                #   extent=(30, 90, vecTbins[ind], vecTbins[ind+1])
                #   )
		#plt.title("Dipole: "+ dipoles[indD]+ "   T = " +str(round(Time_bins[ind],2))+"[h]")
		
#plt.tight_layout()
#fig.suptitle('Ant = ' +sys.argv[2]+', Date = '+sys.argv[1][11:len(sys.argv[1])-14]+', T0 = ' +sys.argv[3]+ ' [h, lst], dT = '+str(dT)+' [h]')
#fig.subplots_adjust(top=0.88)
#if len(sys.argv)==5:
#	plt.show()	
#else:
#	#plt.savefig('plots/Ant' +sys.argv[2]+'_Date'+sys.argv[1][11:len(sys.argv[1])-14]+'_dT'+str(dT)+'.pdf', 		bbox_inches = 'tight',papertype = 'a4')
#	plt.savefig(sys.argv[5])



#print "Time Bins", Time_bins


#print av_data_dict
for i in range(len(vecTbins)-1):
  if len(av_data_dict["Bins"][i].keys()) == 0:
    av_data_dict["Bins"][i]["Not usable"] = ""

hickle.dump(av_data_dict, "dav/averaged_"+sys.argv[1][0:len(sys.argv[1])-4]+ "_ant_" +sys.argv[2]+"_dTh_" +str(dT)+ ".hkl")

#### check the hkl file
#av_data_dict = {}
#def load_hkl(filename):
#	return hickle.load("dav/averaged_"+sys.argv[1][0:len(sys.argv[1])-4]+ "_ant_" +sys.argv[2] +"_dTh_" +str(dT)+ ".hkl")
#avd = load_hkl(sys.argv[1])
#print "New T", avd["Tbins"]
#print "Ant",ant 
#plt.errorbar(avd["Freq"],avd[ant][:], yerr=avd[ant+"_RMS"])
#plt.xlim(42,42.3)
#plt.show()


