import bottleneck as bn
import scipy.signal
import numpy as np
from params import params

def filter(data, size, axis=0):
  # If the input is a masked array, the mask will be lost after filtering

  if not params.median:
    return bn.move_nanmean(data, size, axis=axis)
  else:
    d = np.zeros((data.shape[0], data.shape[1]))
    if size%2 == 0: size -= 1
    if axis == 0:
      for i in range(data.shape[0]):
        d[i] = scipy.signal.medfilt(data[i], size)
    elif axis == 1:
      for i in range(data.shape[1]):
        d[:,i]
        d[:, i] = scipy.signal.medfilt(data[:, i], size)
    else:
      print "Invalid axis", axis, "for filtering"
      exit(1)
    return d

