import sys, hickle

def load_hkl(filename):
  return hickle.load(filename)

data = load_hkl(sys.argv[1])

print "Freqs", "Type", type(data["Freqs"]), "Element Type", type(data["Freqs"][0]), "Shape", data["Freqs"].shape

for i in range(len(data["Bins"])): 

  time_bin = data["Bins"][i]

  print "Bin", i
  for tkey in sorted(time_bin.keys()):
    if tkey == "Time":
     print "LST", time_bin[tkey][0], "->", time_bin[tkey][1], "\tN times", time_bin[tkey][2]
    elif tkey[-3:] == "_SI": print tkey, time_bin[tkey]
    elif tkey == "Not usable": print tkey
    else: print "Key", '"'+tkey+'"', "Type", type(time_bin[tkey]), "Element Type", type(time_bin[tkey][0]), "Shape", time_bin[tkey].shape


