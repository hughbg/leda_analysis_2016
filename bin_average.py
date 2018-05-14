#!/usr/bin/env python

#  Usage: python argv.py outriggers_2018-04-02_07H52M55S.hkl 252 9 0.5 out_hickle_file

import hickle, sys, numpy as np
import matplotlib.pyplot as plt
import math
import bottleneck as bn


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
av_data_dict["Bins"] = []
av_data_dict["Freq"] = d["frequencies"]
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

	data_averaged =  np.zeros((len(vecTbins)-1,len(d["frequencies"])))
	rms_av = np.zeros((len(vecTbins)-1,len(d["frequencies"])))
	

	for ind in range(len(vecTbins)-1):

 		if not use_bins[ind]: continue

		time_bin = {}
		
		indexT = [i for i in range(len(d["lsts"])) if d["lsts"][i]>=vecTbins[ind] and d["lsts"][i]<vecTbins[ind+1]]

		time_bin["Time"] = (d["lsts"][indexT[0]], d["lsts"][indexT[-1]], len(indexT))
		
		data_averaged[ind][:] = np.ma.average(d[ant][indexT][:],axis=0)  	#averaging of the masked array. Masked entries are ignored
 		time_bin[ant] = data_averaged[ind]

		for indnu in range(len(d["frequencies"])):
			#print "indnu", d[ant][indexT,indnu].shape
			flat = bn.move_nanmean(d[ant][indexT,indnu], bp_window_t,axis=0)
			flat = d[ant][indexT,indnu]-flat
			flat = np.ma.filled(flat, np.nan)			# Get rid of masked and NaN, then do std
			flat = flat[np.logical_not(np.isnan(flat))]
			if len(flat) == 0: rms_av[ind][indnu] = np.nan
			else: rms_av[ind][indnu] = np.std(flat)

	#		print "rms",rms_av[ind][indnu]
	#	plt.subplot(len(vecTbins)-1,2,flagsp-2*(len(vecTbins)-1)*indD+indD)
	#	plt.subplot(1,2,flagsp-2*(len(vecTbins)-1)*indD+indD)
		
		time_bin[ant+"_RMS"] = rms_av[ind]
		time_bin[ant+"_weighted"] = data_averaged[ind]/rms_av[ind]**2
		av_data_dict["Bins"].append(time_bin)
	
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


