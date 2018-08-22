# Create s2p files from impedance files

import numpy as np
import cmath

Z_l = 50

files = [ "Zres_252A.txt",  "Zres_252B.txt",  "Zres_254A.txt",  "Zres_254B.txt",  "Zres_255A.txt",  "Zres_255B.txt",  "Zres_256A.txt",  "Zres_256B.txt" ]

for f in files:
  data = np.loadtxt(f)

  f_out = open(f[:-3]+"s2p", "w")
  f_out.write("! Fake S2P created from impedances\n! FREQ.GHZ       S11DB        S11A       S21DB        S21A       S12DB        S12A       S22DB        S22A\n")
  for i in range(21): f_out.write("!\n")

  for d in data:
    freq = d[0]*1e-3		# GHz
    Z_s = complex(d[1], d[3])

    Gamma = (Z_l-Z_s)/(Z_l+Z_s)
    power = 20*cmath.log10(abs(Gamma)).real
    angle = cmath.phase(Gamma)*180/cmath.pi

    f_out.write(" "+( "%8f" %  freq )+"\t"+( "%8f" % power)+"\t"+( "%8f" % angle)+"\t0\t0\t0\t0\t0\t0\n")

  
  f_out.close()

    

  

  


