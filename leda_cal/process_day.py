import ephem, numpy as np
from dpflgr import rfi_flag
from params import params


# RFI flagging and cut out night-time data
# Hold the data and the masks

class ProcessedDay(object):
  def __init__(self, calibrated_data, do_flagging=False, all_lsts=False):

    self.ant_ids = ['252' ] #, '254', '255']
    self.pols = [ "A" ] #, "B" ]

    self.frequencies = calibrated_data["f"]
    self.lsts = calibrated_data["lst"]
    self.utcs = calibrated_data["utc"]

    self.data = {}		# Everything goes in here

    # Report discontinuities in time
    for i in range(1, len(self.lsts)):
        if self.lsts[i]-self.lsts[i-1] > 1/60.0:	# 1 minute
            print "Discontinuity at LST", self.lsts[i], (self.lsts[i]-self.lsts[i-1])*60*60, "seconds"

    self.xlims = (self.frequencies[0], self.frequencies[-1])
    self.ylims = (self.lsts[0], self.lsts[-1])		# Will change if night only
    self.indexes = np.arange(len(self.lsts), dtype=np.int)

    # Get night times, with a border before and after
    self.border_bottom, self.night_bottom, self.night_top, self.border_top = self.calc_night(self.lsts, self.utcs)
    if self.border_bottom: self.night_data_exists = True
    else: self.night_data_exists = False

    if all_lsts:
      for ant in self.ant_ids:
        for p in self.pols:
          print "Processing", ant+p
          T_ant = calibrated_data[ant+p]

          if do_flagging: 
            print "Flagging"
            masks = rfi_flag(T_ant, freqs=self.frequencies)
          else: masks = None

          this_ant = { "data" : T_ant, "masks" : masks }
          self.data[ant+p] = this_ant

    else:
      self.lsts = lst_stamps[self.night_bottom:self.night_top]
      self.utcs = utc_stamps[self.night_bottom:self.night_top]
      self.indexes = indexes[self.night_bottom:self.night_top]
      self.ylims = ( lst_stamps[0], lst_stamps[-1] )

      for ant in self.ant_ids:
        for p in self.pols:
          print "Processing", ant+p
          T_ant = calibrated_data[ant+p]

	  if do_flagging:
            print "Flagging"
            masks = rfi_flag(T_ant[self.border_bottom:self.border_top], freqs=frequencies)
	    masks.chop(self.night_bottom-self.border_bottom, self.night_top-self.border_bottom)		# Knock off border
            T_ant = T_ant[self.night_bottom-self.border_bottom:self.night_top-self.border_bottom]
          else:
            T_ant = T_ant[self.night_bottom: self.night_top]
	    masks = None

          this_ant = { "data" : T_ant, "masks" : masks }
          self.data[ant+p] = this_ant

        # Everything is chopped so that data index 0 is where the night starts. Shift the other pointers.
        self.border_bottom -= self.night_bottom
        self.night_top -= self.night_bottom
        self.border_top -= self.night_bottom
        self.night_bottom = 0

  def report(self):
    print len(self.lsts), "lsts"
    print len(self.utcs), "utcs"
    print "LST range", self.lsts[0], "-", self.lsts[len(self.lsts)-1]
    if self.night_data_exists: print "Night", self.lsts[self.night_bottom], "-", self.lsts[self.night_top-1]
    else: print "Night 0 - 0"

    for ant in self.ant_ids:
      for p in self.pols:
        print ant+p, self.data[ant+p]["data"].shape, type(self.data[ant+p]["data"]), type(self.data[ant+p]["data"][0, 0])
        if self.data[ant+p]["masks"]: self.data[ant+p]["masks"].report()

  def calc_night(self, lsts, utcs):
    if len(lsts) != len(self.utcs):
      raise RuntimeError("LSTs and UTCs are not of the same length")

    ovro_location = ('37.2397808', '-118.2816819', 1183.4839)
    ovro = ephem.Observer(); (ovro.lat, ovro.lon, ovro.elev) = ovro_location
    sun = ephem.Sun()

    gal_center = ephem.FixedBody()  
    gal_center._ra  =  '17 45 40.04'
    gal_center._dec = '-29 00 28.1'
    gal_center.name = "Galactic Center"


    # All these boundary limits are Python style, where the top limit is one above the
    # actual index to use

    night_lsts = []      # night

    # This is the loop that checks every LST to see if the sun/galaxy is up/down
    # and to apply the LST 11:12 restriction.
    for i, d in enumerate(self.utcs):
      ovro.date = d
      sun.compute(ovro)
      gal_center.compute(ovro)


      # This is the big test - all the criteria
      if sun.alt < params.sun_down*np.pi/180 and gal_center.alt < params.galaxy_down*np.pi/180 \
		and params.night_start <= self.lsts[i] and self.lsts[i] <= params.night_end:
        night_lsts.append(i)

    if len(night_lsts) == 0:
      return None, None, None, None
    else: 
      night_lsts.append(night_lsts[-1]+1)		# 1 past the index we want to use

      # Get the border around the time chosen for night, so that we can undo it when flagging done
      window = max(params.st_bp_window_t, params.sc_bp_window_t)

      bottom = max(0, night_lsts[0]-window)
      top = min(len(self.lsts), night_lsts[-1]+window)

      return bottom, night_lsts[0], night_lsts[-1], top		# These are indexes into the data 

