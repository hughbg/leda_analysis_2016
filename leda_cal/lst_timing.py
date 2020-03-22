import numpy as np
import ephem
from params import params


class LST_Timing(object):
  def __init__(self, lsts, utcs):
    self.lsts = lsts
    self.utcs = utcs

  def calc_night(self):
    # Consistency check on different time formats
    # Array of LSTs for the data, and UTCs for the data - each 15s
    # spectrum has these times assigned.
    if len(self.lsts) == 0: return None, None, None, None
    if len(self.lsts) != len(self.utcs):
      raise RuntimeError("LSTs and UTCs are not of the same length")

    # Pyephem setup
    ovro_location = ('37.2397808', '-118.2816819', 1183.4839)
    ovro = ephem.Observer(); (ovro.lat, ovro.lon, ovro.elev) = ovro_location
    sun = ephem.Sun()

    gal_center = ephem.FixedBody()  
    gal_center._ra  =  '17 45 40.04'
    gal_center._dec = '-29 00 28.1'
    gal_center.name = "Galactic Center"


    night_lsts = []      # An array of indexes into the data where is night.
			# The same indexes apply to LST and UTC arrays.

    # This is the loop that checks every LST to see if the sun/galaxy is up/down
    # and to apply the LST hard restriction. The index is "i", which is saved if night.
    for i, d in enumerate(self.utcs):
      # Find where the sun and galaxy are at this time.
      ovro.date = d
      sun.compute(ovro)
      gal_center.compute(ovro)


      # This is the big test - all the criteria. Check the altitude of sun/galaxy
      # and if LST is out of range. Limits are set in params.py.
      # Problem if night is not contiguous chunk - could it not be?
      if sun.alt < params.sun_down*np.pi/180 and gal_center.alt < params.galaxy_down*np.pi/180 \
		and params.night_start <= self.lsts[i] and self.lsts[i] <= params.night_end:
        night_lsts.append(i)

    if len(night_lsts) == 0:	# Means there is no night time in the data.
      return None, None, None, None
    else: 

      # Get the border around the time chosen for night, so that we can added it when flagging done.
      # The border is included when flagging is done, but stripped out when the data is saved.
      window = max(params.st_bp_window_t, params.sc_bp_window_t)

      # bottom/top is where the border is that wraps the night-time.
      # max/min used so that we don't go outside the array.
      bottom = max(0, night_lsts[0]-window)
      top = min(len(self.lsts)-1, night_lsts[-1]+window)

      return bottom, night_lsts[0], night_lsts[-1], top		# These are indexes into the data 
   
  def align(self, data):   # This is not working. Meant to allign the data to 15s because sometimes it's off.
    num_lsts_in_day = (24*60*60)/15
    thirteen_seconds = 13.0/60/60		# As decimal number 0...24
    fifteen_seconds = 15.0/60/60
    seventeen_seconds = 17.0/60/60

    #print 60*60*(lsts[-1]-lsts[0])/(len(lsts)-1); exit()

    new_data = np.zeros((num_lsts_in_day, data.shape[1]))
    new_lsts = np.zeros(num_lsts_in_day)
    new_utcs = np.zeros(num_lsts_in_day)
    counts = np.zeros(num_lsts_in_day)

    # Step through the times and bin data
    for i in range(len(self.lsts)):
      bin_index = self.lst2bin(self.lsts[i])
      new_data[bin_index, :] += data[i, :]
      counts[bin_index] += 1

      new_lsts[bin_index] = self.bin2lst(bin_index)

    for i in range(num_lsts_in_day):
      if counts[i] > 1: new_data[i, :] /= counts[i]

    return new_lsts, new_data

  def lst2bin(self, lst):
    return int(lst*60*60/15)

  def bin2lst(self, bin_num):
    return float(bin_num)*15/60/60


