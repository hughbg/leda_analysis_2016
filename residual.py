# Probably should use 05 script or merge them.

import settings
import numpy as np
import matplotlib.pyplot as plt

spectrum = np.loadtxt(settings.ledaSpec)
frequencies = np.loadtxt(settings.ledaFreqs)

spectrum1 = spectrum[spectrum!=0.0]
frequencies = frequencies[spectrum!=0.0]
spectrum = spectrum1


p = [ 1.0, 2.0, 3.0, 4.0 ]
for line in open(settings.outdir+"/1-stats.dat"):
    l = line.split()
    if len(l) > 0:
        if l[0] in [ "4", "5", "6", "7" ]:
            p[int(l[0])-4] = float(l[1])

foreground = np.zeros(len(frequencies))

f0 = 60

for j in range(len(frequencies)):
    f = frequencies[j]
    for n in range(len(p)):
        foreground[j] += p[n]*np.power(np.log10(f/f0), n)
    foreground[j] = np.power(10, foreground[j])

plt.ylabel("Temp")
plt.xlabel("Frequency [MHz]")
plt.plot(frequencies, spectrum, label="Spectrum")
plt.plot(frequencies, foreground, label="Foreground")
plt.title("Spectrum and Foreground")
plt.legend()
plt.show()

    
     
