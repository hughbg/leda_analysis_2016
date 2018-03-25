import hickle, sys, numpy as np
import matplotlib.pyplot as plt

def load_hkl(filename):
  return hickle.load(filename)

def report_hkl(data):
  print "Items present:"
  for key in sorted(data.keys()): 
    try:
      if len(data[key].shape) == 1: print "Key", '"'+key+'"', "Type", type(data[key]), "Element Type", type(data[key][0]), "Shape", data[key].shape
      else: print "Key", '"'+key+'"', "Type", type(data[key]), "Element Type", type(data[key][0][0]), "Shape", data[key].shape
    except: pass
  print "\nStart time", data["utcs"][0]
  print "Options", data["options"]

  # Sometimes you have to relpace masked data with a value
  plt.imshow(np.ma.filled(data["252A"], 0), clim=(1000,10000), aspect="auto")
  plt.show()

d = load_hkl(sys.argv[1])
report_hkl(d)
