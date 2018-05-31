class Params(object):
  def __init__(self):
    """
    # Default params, for restoring 
    self.sun_down = -15		# In degrees. Below this they are down.
    self.galaxy_down = -15
    self.do_sum_threshold = True
    self.do_sigma_clip = True
    self.do_dtv_flagging = True
    self.thr_f = 0.2
    self.thr_t = 0.2
    self.rho = 1.5
    self.st_bp_window_f = 16			# Sliding windows in time and frequency
    self.st_bp_window_t = 16
    self.max_frac_f = 0.5
    self.max_frac_t = 0.5
    self.thr_max = 10
    self.thr_min = 0
    self.scales = [1, 2, 4, 8, 16, 32, 64]
    self.median = False			# Use scipy median filter instead of bottleneck mean for sum threshold
    self.sc_bp_window_f = self.st_bp_window_f			# Sliding windows in time and frequency
    self.sc_bp_window_t = self.st_bp_window_t
    self.sigma = 4						# Clip above/below this
    self.un_bp_window_t = self.st_bp_window_t		# Window in time for flattening 1 channel to get RMS of channel
    self.ds9_clip = 4		
    self.stats_bp_window_f = self.st_bp_window_f 		# For flattening
    self.stats_bp_window_t = self.st_bp_window_t
    self.histogram_length = 5000			# Length for gauss fitting histogram of waterfall plot data

    self.dtv_frequencies = [ 54.31, 60.31, 66.31, 76.31, 82.31 ]	# Pilot tone
    """

    # For day/night selection
    self.sun_down = -15		# In degrees. Below this they are down.
    self.galaxy_down = -15

    # Select flagging algorithms
    self.do_sum_threshold = True
    self.do_sigma_clip = True
    self.do_dtv_flagging = True

    # Sum threshold flagging params. Some description in the code describes these.
    self.thr_f = 0.2
    self.thr_t = 0.2
    self.rho = 1.5
    self.st_bp_window_f = 16			# Sliding windows in time and frequency
    self.st_bp_window_t = 16
    self.max_frac_f = 0.5
    self.max_frac_t = 0.5
    self.thr_max = 10
    self.thr_min = 0
    self.scales = [1, 2, 4, 8, 16, 32, 64]
    self.median = False			# Use scipy median filter instead of bottleneck mean for sum threshold

    # Parameters for flagging using sigma clipping
    self.sc_bp_window_f = self.st_bp_window_f			# Sliding windows in time and frequency
    self.sc_bp_window_t = self.st_bp_window_t
    self.sc_passes = 3					# Number of passes to make through the clipper
    self.sigma = 3						# Clip above/below this

    # Uncertainties. 
    self.un_bp_window_t = self.st_bp_window_t		# Window in time for flattening 1 channel to get RMS of channel

    # DS9 clip. If flagging is not being done, it is necessary to clip the data for visualizing in DS9. 
    # There are very large peaks in the data, and if these are included in the FITS file, the huge range of the data
    # makes it difficult to see any detail outside of the peaks. The clip is a sigma clip _without_ flattening first.
    self.ds9_clip = 4		


    # Parameters for statistics gathering. Windows for statistics are independent of flagging windows because
    # we need to gather statistics from data that is not flagged. Also we may want to vary the flagging independently
    # of the statistics.	
    self.stats_bp_window_f = self.st_bp_window_f 		# For flattening
    self.stats_bp_window_t = self.st_bp_window_t
    self.histogram_length = 5000			# Length for gauss fitting histogram of waterfall plot data

    self.dtv_frequencies = [ 54.31, 60.31, 66.31, 76.31, 82.31 ]	# Pilot tone
    
  def __repr__(self):
    order = ['sun_down', 'galaxy_down', 'do_sum_threshold', 'do_sigma_clip', 'do_dtv_flagging', 
             'thr_f', 'thr_t', 'rho', 'st_bp_window_f', 'st_bp_window_t', 'max_frac_f', 
             'max_frac_t', 'thr_max', 'thr_min', 'scales', 'median', 'sc_bp_window_f', 
             'sc_bp_window_t', 'sigma', 'un_bp_window_t', 'ds9_clip', 
             'stats_bp_window_f', 'stats_bp_window_t', 'histogram_length', 'dtv_frequencies']
    
    output = ""
    for key in order:
      output += "%s: %s\n" % (key, getattr(self, key, None))
    return output


params = Params()
