import numpy as np

data = np.loadtxt("leda.lna.s11.cable.de-embedded.csv", delimiter=",", skiprows=1)

# Create fake s2p file
def s2p(ant, col_dB, col_deg):
  f = open("leda.lna.s11.cable.de-embedded"+"_"+ant+".s2p", "w")
  f.write("! Fake S2P LNA data\n! FREQ.GHZ       S11DB        S11A       S21DB        S21A       S12DB        S12A       S22DB        S22A\n")

  for d in data:
    freq = d[0]*1e-3		# GHz
    f.write(" "+( "%8f" %  freq )+"\t"+( "%8f" % d[col_dB])+"\t"+( "%8f" % d[col_deg])+"\t0\t0\t0\t0\t0\t0\n")

  f.close()



s2p("252A", 1, 2)
s2p("252B", 3, 4)
s2p("254A", 5, 6)
s2p("254B", 7, 8)
s2p("255A", 9, 10)
s2p("255B", 11, 12)
s2p("256B", 13, 14)

