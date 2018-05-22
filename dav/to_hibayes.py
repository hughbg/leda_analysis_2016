import sys, hickle, numpy as np


if len(sys.argv) != 4:
  print "Usage:", sys.argv[0], " hkl_file  which_ant which_bin"
  print "	hkl_file:	from bin_average.py"
  print "	which_ant:	like 252A"
  print "	which_bin:	one of the numbered time bins in the file"
  exit(1)

data =  hickle.load(sys.argv[1])

np.savetxt("frequencies.dat", data["Freqs"][19:])

if int(sys.argv[3]) >= len(data["Bins"]):
  print "Bin does not exist"
  exit(1)

time_bin = data["Bins"][int(sys.argv[3])]

for tkey in sorted(time_bin.keys()):
  if tkey == "Not usable": print tkey
  if tkey == sys.argv[2]: np.savetxt("spectrum.dat", time_bin[tkey][19:])
  if tkey == sys.argv[2]+"_RMS": np.savetxt("spectrum_errs.dat", time_bin[tkey][19:])


