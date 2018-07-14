import numpy as np
import cmath
import scipy.interpolate


# Constants
Z_0 = 50
Y_0 = 1/50.0
T_amb_value = 300

def degree2radian(a):
  return a*cmath.pi/180

# T_hot and T_ambient are not loaded from files.
def generate_T_amb_hot(length):
  T_amb = np.full(length, T_amb_value)		# This is an array, so then T_hot will be an array
  Q = 14.9-6.95
  T_hot = T_amb*(10**(Q/10)+1)

  return T_amb, T_hot


# Load data from an s2p file and convert columns into Gamma values and interpolate the values onto LEDA frequencies
def load_Gamma(fname, leda_frequencies):
  col_S22DB = 7			# Columns where S22DB and S22A are
  col_S22A = 8

  f = open(fname)
  data = np.loadtxt(f.readlines()[23:])		# Load data and ignore header
  f.close()

  values = [ cmath.rect(10**(row[col_S22DB]/20), degree2radian(row[col_S22A])) for row in data ]	# Convert to complex
  s2p_frequencies = data[:, 0]*1e9			# Work in Hz, convert to Hz

  # Now got the values from the file, but need to to interpolate values for LEDA channels.
  # Generate the function for interpolation
  function = scipy.interpolate.interp1d(s2p_frequencies, values)

  # Now generate values for LEDA channels. However can't interpolate outside the range of the s2p frequencies
  # so must check.

  data = np.zeros(len(leda_frequencies), dtype=np.complex64)
  for i in range(len(leda_frequencies)):
    leda_frequency = leda_frequencies[i]*1e6
    if leda_frequency < s2p_frequencies[0] or s2p_frequencies[-1] < leda_frequency: data[i] = complex(0,0)
    else: data[i] = function(leda_frequencies[i]*1e6)

  return data


# This function does the work (for a single frequency). Most values are scalars except for P_meas which
# is a list with 4 measurements

def do_calculation(P_hot, P_cold, T_hot, T_amb, Gamma_hot, Gamma_cold, Gamma_lna, P_meas, freq):
  if len(P_meas) != 4:
    print "Error, expecting 4 measurements OFF, SHORT, 47pf, 66pf",
    exit(1)
  T_hot = 2185
  print "P_hot", P_hot, "P_cold", P_cold, "T_hot", T_hot, "T_amb", T_amb, "Gamma_hot", Gamma_hot, "Gamma_cold", Gamma_cold, "Gamma_lna", Gamma_lna, "P_meas", P_meas
  

  # Firstly sort out Gamma_s for the 4 measurements

  Z_47 = complex(0, -1/(2*cmath.pi*freq*1e6*47e-12))	# Z = -j/(angular frequency * C)
  Gamma_s_47 = (Z_47-Z_0)/(Z_47+Z_0)
  Z_66 = complex(0, -1/(2*cmath.pi*freq*1e6*66e-12))
  Gamma_s_66 = (Z_66-Z_0)/(Z_66+Z_0)
  Gamma_s = [ 0, -1, 1, Gamma_s_47  ]		# OFF, OPEN, SHORT, 47

  # Start calculating

  Gamma_ns = (Gamma_hot+Gamma_cold)/2
  S_P_T = (P_hot-P_cold)/(T_hot-T_amb) * abs(1-Gamma_lna*Gamma_ns)**2/(1-abs(Gamma_ns)**2)
  print "S_P_T", S_P_T

  # Should be calculated as: Gamma_s = Gamma_T* (S21**2)
  # where Gamma _T is the termination reflection coefficient.
  # These values work for the test files.
  Gamma_s =   [0, complex(0.836,-0.38), complex(-0.836, 0.38), complex(-0.406, -0.821) ]
  print "Gamma_s", Gamma_s

  # Build T, X, matrices

  T = []
  X = []
  for i in range(len(P_meas)):
    T.append([ ((P_meas[i]/S_P_T)*abs(1-Gamma_lna*Gamma_s[i])**2 - T_amb*(1-abs(Gamma_s[i])**2)) / T_amb ])	# T_R_i
  
    x = abs(1-Gamma_lna*Gamma_s[i])**2
    y = 1-abs(Gamma_s[i])**2

    #print (T[0][0]+y)*T_amb/P_meas[i], x/S_P_T; exit()

    X.append([ (1-abs(Gamma_s[i])**2), abs(1-Gamma_s[i])**2, abs(1+Gamma_s[i])**2, 2*Gamma_s[i].imag ] )	# X_row_i
  
  print "T",T
  print "X", X
  T = np.matrix(T)
  X = np.matrix(X)

  # Get vector C, using matrix ops https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.matrix.html

  #C = ((X.H*X).I)*(X.H*T)

  C = X.I*T				# we have square matrix so do it the easy way

  print C
  
  # Extract result values
  R_N = C[1]/Y_0
  B_opt = Y_0*C[3]/C[1]					# Note vectors are indexed base 0 in Python
  G_opt = Y_0/C[1]*cmath.sqrt(C[1]*C[2]-C[3]**2)
  T_min = T_amb*(C[0]+2*cmath.sqrt(C[1]*C[2]-C[3]**2))

  return ( R_N, B_opt, G_opt, T_min )


############### MAIN PROGRAM
# Deals mostly with loading the data


# Load the files into arrays. Get a list of frequencies first, they are in MHz. 
frequencies = np.loadtxt("ant_254A.LNA0.yf346-7.off.1.2018-05-24_14-29-09.dat")[:, 0]   # Better be the same frequencies in all the P_ files!
										        # Keep a frequency list for later on.
array_P_hot = np.loadtxt("ant_254A.LNA0.yf346-7.on.1.2018-05-24_14-30-06.dat")[:, 1]	# The [:, 1] term throws away the frequency column
array_P_cold = np.loadtxt("ant_254A.LNA0.yf346-7.off.1.2018-05-24_14-29-09.dat")[:, 1]
array_T_amb, array_T_hot = generate_T_amb_hot(len(frequencies))				# All constants
array_Gamma_hot = load_Gamma("346-7bw3.off.s11.s2p", frequencies)
array_Gamma_cold = load_Gamma("346-7bw3.on.s11.s2p", frequencies)
array_Gamma_lna = load_Gamma("254a.lna.rl.18may13.s2p", frequencies)
	# Gammas_ns is calculated later

# Load P_meas, there will be 4 measurements in 4 files placed in the list below.
# The order MUST be OFF, OPEN,  SHORT, 47 to match Gamma_s above.
files = [ "ant_254A.LNA0.yf346-7.off.1.2018-05-24_14-29-09.dat",  "ant_254A.LNA0.2p0m.OPEN.1.2018-05-24_14-33-49.dat", "ant_254A.LNA0.2p0m.SHORT.1.2018-05-24_14-35-44.dat",  "ant_254A.LNA0.2p0m.47pf.1.2018-05-24_14-37-37.dat" ]

array_P_meas = [ None for i in range(len(files)) ]		# Gather the P_meas data into a list
for i in range(len(files)):
  array_P_meas[i] = np.loadtxt(files[i])[:, 1]


# Sanity check - check lengths all the same
all_data = [ array_P_hot, array_P_cold, array_T_hot, array_T_amb, array_Gamma_hot, array_Gamma_cold, array_Gamma_lna ]
for i in range(1, len(all_data)):
  if len(all_data[i]) != len(all_data[0]):
    print "Error: input arrays are not all the same length, got",  len(all_data[i]), len(all_data[0])
    exit(1)
for i in range(len(array_P_meas)):
  if len(array_P_meas[i]) != len(all_data[0]):
    print "Error: input P_meas arrays are not the right length"
    exit(1)



array_P_meas = np.array(array_P_meas)	

# Loop over frequencies from 30-88MHz
for i in [ 2500 ]: #range(1250, 3666):
  print "Freq", frequencies[i],
  try:
     R_N, B_opt, G_opt, T_min  = do_calculation(array_P_hot[i], array_P_cold[i], array_T_hot[i], array_T_amb[i], array_Gamma_hot[i], 
  	array_Gamma_cold[i], array_Gamma_lna[i], array_P_meas[:, i], frequencies[i])
     print "R_N", R_N, "B_opt", B_opt, "G_opt", G_opt, "T_min", T_min
  except:
    print "failed"


