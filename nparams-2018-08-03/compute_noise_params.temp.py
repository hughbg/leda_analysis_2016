#!/usr/bin/env python
# coding: utf-8
"""
# compute_nparams.py

Compute noise parameters for a device-under-test using Edward's method.
"""

import numpy as np
import pylab as plt
from scipy.interpolate import interp1d as interp

# Constants
Z_0 = 50.0
Y_0 = 1/50.0
T_0 = 290.0

def to_complex(data, linear=True, radians=False):
    """ Convert amp / phase data to complex."""
    r = data[:, 0]
    q = data[:, 1]

    if linear is False:
        r = 10**(r / 20)
    if radians is False:
        q = np.deg2rad(q)

    x = r * np.cos(q)
    y = r * np.sin(q)
    return x + 1j* y

class Spectrum(object):
    """ Class for reading spectrometer dumps

    Provides interpolation function and basic plotting of data
    """
    def __init__(self, filename):
        d = np.loadtxt(filename)
        self.filename = filename
        self.f = d[:, 0]
        self.d = interp(d[:,0], d[:, 1])
    def plot(self, log=True):
        if log:
            plt.plot(self.f, 10*np.log10(self.d(self.f)))
        else:
            plt.plot(self.f, self.d(self.f))
        plt.xlabel("Frequency [MHz]")

def read_spectrum(filename):
    return Spectrum(filename)

class Sparam(object):
    """ Basic class for interpolating and plotting S-parameter data """

    def __init__(self, f=None, s11=None, s21=None, s12=None, s22=None):
        self.f    = f
        self.s11  = s11
        self.s21  = s21
        self.s22  = s22
        self.s12  = s12

    def plot_s11(self):
        """ Plot S11 data """
        s11 = self.s11(self.f)
        s11_mag = 20*np.log10(np.abs(s11))
        s11_phs = np.rad2deg(np.angle(s11))
        plt.subplot(1,2,1)
        plt.plot(self.f, s11_mag)
        plt.subplot(1,2,2)
        plt.plot(self.f, s11_phs)

    def plot_s21(self):
        """ Plot S21 data """
        s21 = self.s21(self.f)
        s21_mag = 20*np.log10(np.abs(s21))
        s21_phs = np.rad2deg(np.angle(s21))
        plt.subplot(1,2,1)
        plt.plot(self.f, s21_mag)
        plt.subplot(1,2,2)
        plt.plot(self.f, s21_phs)


def read_anritsu_s11(filename):
    """ Read data from Anristu VNA and return instance of Sparam class """
    d = np.genfromtxt(filename, delimiter=',', skip_header=8, skip_footer=1)
    f_mhz = d[:, 0] / 1e6
    s11 = interp(f_mhz, to_complex(d[:,1:], linear=False, radians=False))
    S = Sparam(f_mhz, s11)
    return S

def read_cable_sparams(filename):
    """ Read data from cable CSV files and return instance of Sparam class """
    d = np.genfromtxt(filename, delimiter=',', skip_header=1, skip_footer=0)
    f = d[:, 0]
    s11 = interp(f, to_complex(d[:,1:3], linear=False, radians=False))
    s21 = interp(f, to_complex(d[:,3:5], linear=False, radians=False))
    s22 = interp(f, to_complex(d[:,5:7], linear=False, radians=False))
    S = Sparam(f, s11=s11, s21=s21, s22=s22)
    return S

def read_uu_sparams(filename):
    """ Read data from cable CSV files and return instance of Sparam class """
    d = np.genfromtxt(filename, skip_header=0, skip_footer=0)
    print d
    f = d[:, 0]

    s21 = interp(f, to_complex(d[:,1:3], linear=False, radians=False))
    S = Sparam(f, s21=s21)
    return S

def read_s2p_s11(filename, s11_col=7):
    """ Read data from S2P VNA file and return instance of Sparam class.

    Only reads S11 data, defaults to column 7 (sometimes it is column 1)
    """
    c = s11_col
    with open(filename) as fh:
        data = np.loadtxt(fh.readlines()[23:])		# Load data and ignore header
    f_mhz  = data[:, 0] * 1e3 # To MHz
    s11    = to_complex(data[:, c:c+2], linear=False, radians=False)
    s11_i  = interp(f_mhz, s11)
    return Sparam(f_mhz, s11_i)

def generate_T_amb_hot(length):
    """ T_hot and T_ambient are not loaded from files. """
    T_amb = np.full(length, 32+273.15)	# This is an array, so then T_hot will be an array
    Q = 14.9-6.95
    T_hot = T_amb*(10**(Q/10)+1)

    return T_amb, T_hot

def compute_noise_coeffs(f_mhz, P_hot, P_cold, T_hot, T_cold, GM_hot, GM_cold, GM_lna, GM_meas, P_meas):
    """ Apply Edward's method to compute noise vector C for a device-under-test

    Args:
        f_mhz (np.array):   Array of frequency values, in MHz.
        P_hot (np.array):   Measured power for hot reference load,
                            in counts (uncalibrated)
        P_cold (np.array):  Measure power for ambient reference load,
                            in counts (uncalibrated)
        T_hot (np.array):   Noise temperature for hot reference load, in K.
        T_cold (np.array):  Noise temperature for ambient reference load, in K.
        GM_hot (np.array):  Reflection coefficient of hot load (complex linear)
        GM_cold (np.array): Reflection coefficient of cold load (complex linear)
        GM_lna (np.array):  Reflection coefficient of LNA (or device-under-test)
        GM_meas (list of np.arrays): List of refl. coeffs. for all references
        P_meas  (list of np.arrays): List of uncalibrated power measurements
                                     for all reference sources.

    Returns:
        Tuple of noise vectors, [C0. C1, C2, C3]

    Notes:
        Returns noise vector [C0. C1, C2, C3]; use C_to_nparams()
        to convert to standard noise parameters.

    """
    # For normalization, use average of hot and ambient refl coefficients
    GM_ns = (GM_hot+GM_cold)/2

    # Compute S_P_T, which is used to normalize measured power to temperature
    S_P_T = (P_hot-P_cold)/(T_hot-T_cold) * np.abs(1-GM_lna*GM_ns)**2/(1-np.abs(GM_ns)**2)

    # Generate a 4xN array (N is frequency axis which we loop over)
    X0 = 1.0-np.abs(GM_meas)**2
    X1 = np.abs(1.0-GM_meas)**2
    X2 = np.abs(1.0+GM_meas)**2
    X3 = 4*GM_meas.imag
    X = np.hstack((np.expand_dims(X0, 1),
                   np.expand_dims(X1, 1),
                   np.expand_dims(X2, 1),
                   np.expand_dims(X3, 1)))

    T_Ri = ((P_meas / S_P_T) * abs(1-GM_lna*GM_meas)**2) / T_cold

    # Now compute C params, looping over frequency axis
    C = []
    for ii in range(T_Ri.shape[1]):
        Tm = np.matrix(T_Ri[..., ii]).T
        Xm = np.matrix(X[..., ii])
        #print Tm.shape, Xm.shape
        #Cm = Xm.I*Tm # we have square matrix so do it the easy way
        Cm = ((Xm.H*Xm).I)*(Xm.H*Tm)
        C.append(Cm)
    return np.array(C).squeeze()

def C_to_nparams(C):
    """ Convert C-vector noise quantities to regular noise params """
    # Note vectors are indexed base 0 in Python
    R_N = C[:, 1]/Y_0
    B_opt = Y_0*C[:, 3]/C[:, 1]
    G_opt = Y_0/C[:, 1]*np.sqrt(C[:, 1]*C[:, 2]-C[:, 3]**2)
    F_min = C[:, 0] + 2*np.sqrt(C[:, 1]*C[:, 2]-C[:, 3]**2)
    return [ R_N, B_opt, G_opt, F_min ]

def y_factor(T_hot, T_cold, P_hot, P_cold):
    """ Compute Y-factor via standard method

    Args:
        T_hot: Hot load temperature in K
        T_cold: Cold load temperature in K
        P_hot: Power measured from HOT, linear units
        P_cold: Power measured from COLD, linear units
    """
    y = P_hot / P_cold
    T_dut = ((T_hot) - y * T_cold) / (y - 1)
    return T_dut

def NT_to_NF(T):
    """ Convert noise temperature to noise figure """
    return 1 + T / T_0

def NF_to_NT(F):
    """ Convert noise figure to noise temperature """
    return T_0 * (F - 1)

def Y_to_GM(Y):
    """ Convert complex admittance Y to reflection coeff Gamma """
    GM = (1.0 - Y * Z_0) / (1.0 + Y * Z_0)
    return GM

def compute_NF(F_min, R_N, GM_S, GM_opt):
    """ Compute noise figure F given F_min, R_N, GM_S, and GM_opt """
    a0 = (4.0 * R_N / Z_0) * (np.abs(GM_S - GM_opt)**2)
    a1 = (1.0 - np.abs(GM_S)**2) * np.abs(1.0 + GM_opt)**2
    F = F_min + a0 * a1
    return F

def compute_T_nw(F_min, R_N, B_opt, G_opt, GM_S):
    """ Compute noisewave for given """
    Y_opt  = G_opt + 1j*B_opt
    GM_opt = Y_to_GM(Y_opt)
    F      = compute_NF(F_min, R_N, GM_S, GM_opt)
    T      =  NF_to_NT(F)
    return T

def compute_G_S(S11_S, S11_dut):
    """ Compute Transducer gain quantity G_S """
    G_S = (1.0 - np.abs(S11_S)) / np.abs(1.0 - S11_S * S11_dut)
    return G_S

def apply_calibration(P_3ss, T_fe_hot, T_fe_cold,
                      T_rx_cal, T_rx_ant,
                      S11_ant, S11_lna):
    """ Apply calibration formalism to 3SS data

    Args:
        P_3ss (np.array): (P_sky - P_cold) / (P_hot - P_cold)
        T_fe_hot (np.array): Frontend noise diode hot reference
        T_fe_cold (np.array): Frontend noise diode cold reference
        T_rx_cal (np.array): Receiver temperature when connected to cal source
        T_rx_ant (np.array): Receiver temperature when connected to antenna
        S11_ant (np.array): Reflection coefficient of antenna
        S11_lna (np.array): Reflection coefficient of LNA
    """
    G_S = compute_G_S(S11_ant, S11_lna)
    A = (T_fe_hot - T_fe_cold)
    B = T_fe_cold + T_rx_cal -  T_rx_ant * G_S

    #plt.subplot(411)
    #plt.plot(T_fe_cold)
    #plt.subplot(412)
    #plt.plot(T_fe_hot)
    #plt.subplot(413)
    #plt.plot(T_rx_cal)
    #plt.subplot(414)
    #plt.plot(T_rx_ant)
    #plt.show()
    #
    #plt.subplot(311)
    #plt.plot(A)
    #plt.subplot(312)
    #plt.plot(B)
    #plt.subplot(313)
    #plt.plot(1.0 / G_S)
    #plt.show()

    T_sky = (1.0 / G_S) * (A * P_3ss + B)
    return T_sky

if __name__ == "__main__":


    # Set frequency ranges to interpolate data onto
    f_mhz = np.linspace(30, 87.5, 201)
    testpoint = 'LNA0'

    ###################################
    ##   READ VNA AND SPECTRUM DATA  ##
    ###################################

## Two sets of noise parameters for each signal path.  One measured at input to switch.  One at input to LNA.
## Measurements explicitly tied to noise parameter measurements:
##  (i)   2m and 0.9m cables w/ 5 different end loads:  open, short, capacitance 1, capacitance 2, and 50 ohm termination
##  (ii)  Injected power by Keysight noise source, two states, ON and OFF
##
## For SW0 test point measurements are made in a 5-5-5 sec cycle.  Only one of the 5 sec intervals delivers data for the
## 0.9 and 2m cables.  The other two report the power of the hot and the cold diode paths.  The latter data can be used to cancel
## LNA variability.  NOTE:  The system varies slowly.  It SHOULD not be necessary to use EVERY 5 sec record for the hot and cold
## diode states.
##
## For LNA0 test point, each measurement generates three files, one for each 5 sec.  These can be averaged together for use.


    if testpoint == 'SW0':
        # Load lab measurements of open, short, termination and capacitor reflection coefficients   ???
        # LJG NOTE: Need to examine origin of cable S11s and to update w/ cascaded btw-series adapters.
        c47 = read_anritsu_s11('cal_data/c47pF_0702.csv')
        c66 = read_anritsu_s11('cal_data/c66pF_0702.csv')
        o   = read_anritsu_s11('cal_data/open0702.csv')
        s   = read_anritsu_s11('cal_data/short0702.csv')
        l   = read_anritsu_s11('cal_data/R50p9_0702.csv')

        # Load 0.9 and 2m cable S-params
        cable_0p9m = read_cable_sparams('cal_data/cable_0p9.csv')
        cable_2m   = read_cable_sparams('cal_data/cable_2m.csv')

        # Load Keysight 346B noise source S11 values
        s2p_hot  = read_s2p_s11('cal_data/346-7bw3.on.s11.s2p')
        s2p_cold = read_s2p_s11('cal_data/346-7bw3.off.s11.s2p')

        # Load LNA S11 values
        s2p_lna  = read_s2p_s11('cal_data/254A/254a.lna.rl.18may13.s2p')
    
        # Read balun S11 measurements
        s2p_ant = read_s2p_s11('cal_data/254A/4a.m.1.s2p', s11_col=1)

        # Now load uncalibrated spectra corresponding to reference sources
        P_2m_open    = read_spectrum('cal_data/254A/ant_254A.SW0.2p0m.OPEN.skypath.dat')
        P_2m_short   = read_spectrum('cal_data/254A/ant_254A.SW0.2p0m.SHORT.skypath.dat')
        P_2m_c47     = read_spectrum('cal_data/254A/ant_254A.SW0.2p0m.47pf.skypath.dat')
        P_2m_c66     = read_spectrum('cal_data/254A/ant_254A.SW0.2p0m.66pf.skypath.dat')
        P_0p9m_open  = read_spectrum('cal_data/254A/ant_254A.SW0.0p9m.OPEN.skypath.dat')
        P_0p9m_short = read_spectrum('cal_data/254A/ant_254A.SW0.0p9m.SHORT.skypath.dat')
        P_0p9m_c47   = read_spectrum('cal_data/254A/ant_254A.SW0.0p9m.47pf.skypath.dat')
        P_0p9m_c66   = read_spectrum('cal_data/254A/ant_254A.SW0.0p9m.66pf.skypath.dat')

        # Load / compute spectra and temperature for hot and ambient reference sources
        T_cold, T_hot = generate_T_amb_hot(len(f_mhz))                                        ## ???
        P_hot      = read_spectrum('cal_data/254A/ant_254A.SW0.yf346-7.on.skypath.dat')
        P_cold     = read_spectrum('cal_data/254A/ant_254A.SW0.yf346-7.off.skypath.dat')


        # Load noise diode states - select hot and cold from one of many configurations.
        P_fe_cold    = read_spectrum('cal_data/254A/ant_254A.SW0.0p9m.OPEN.coldpath.dat')
        P_fe_hot     = read_spectrum('cal_data/254A/ant_254A.SW0.0p9m.OPEN.hotpath.dat')

    elif testpoint == 'LNA0':
            # Load lab measurements of open, short, termination and capacitor reflection coefficients
        
            c47 = read_anritsu_s11('cal_data/c47pF_0702.csv')
            c66 = read_anritsu_s11('cal_data/c66pF_0702.csv')
            o   = read_anritsu_s11('cal_data/open0702.csv')
            s   = read_anritsu_s11('cal_data/short0702.csv')
            l   = read_anritsu_s11('cal_data/R50p9_0702.csv')

            # Load 0.9 and 2m cable S-params
            cable_0p9m = read_cable_sparams('cal_data/cable_0p9.csv')
            cable_2m   = read_cable_sparams('cal_data/cable_2m.csv')

            # Load Keysight 346B noise source S11 values
            s2p_hot  = read_s2p_s11('cal_data/346-7bw3.on.s11.s2p')
            s2p_cold = read_s2p_s11('cal_data/346-7bw3.off.s11.s2p')

            # Load LNA S11 values
            s2p_lna  = read_s2p_s11('cal_data/254A/254a.lna.rl.18may13.s2p')

            # Read balun S11 measurements
            s2p_ant = read_s2p_s11('cal_data/254A/4a.m.1.s2p', s11_col=1)

            # Now load uncalibrated spectra corresponding to reference sources
            P_2m_open    = read_spectrum('cal_data/254A.LNA0/ant_254A.LNA0.2p0m.OPEN.dat')
            P_2m_short   = read_spectrum('cal_data/254A.LNA0/ant_254A.LNA0.2p0m.SHORT.dat')
            P_2m_c47     = read_spectrum('cal_data/254A.LNA0/ant_254A.LNA0.2p0m.47pf.dat')
            P_2m_c66     = read_spectrum('cal_data/254A.LNA0/ant_254A.LNA0.2p0m.66pf.dat')
            P_0p9m_open  = read_spectrum('cal_data/254A.LNA0/ant_254A.LNA0.0p9m.OPEN.dat')
            P_0p9m_short = read_spectrum('cal_data/254A.LNA0/ant_254A.LNA0.0p9m.SHORT.dat')
            P_0p9m_c47   = read_spectrum('cal_data/254A.LNA0/ant_254A.LNA0.0p9m.47pf.dat')
            P_0p9m_c66   = read_spectrum('cal_data/254A.LNA0/ant_254A.LNA0.0p9m.66pf.dat')

            # Load / compute spectra and temperature for hot and ambient reference sources
            T_cold, T_hot = generate_T_amb_hot(len(f_mhz))                                    ##  ???
            P_hot      = read_spectrum('cal_data/254A.LNA0/ant_254A.LNA0.yf346-7.on.dat')
            P_cold     = read_spectrum('cal_data/254A.LNA0/ant_254A.LNA0.yf346-7.off.dat')

            # Load noise diode states - select hot and cold from one of many configurations.   ??? ERROR ???
            P_fe_cold    = read_spectrum('cal_data/254A.LNA0/ant_254A.coldpath.dat')
            P_fe_hot     = read_spectrum('cal_data/254A.LNA0/ant_254A.hotpath.dat')

    # Compute reflection coefficient for cable + O/S/L standards
    cable_2m_open    = o.s11(f_mhz) * cable_2m.s21(f_mhz)**2
    cable_2m_short   = s.s11(f_mhz) * cable_2m.s21(f_mhz)**2
    cable_2m_load    = l.s11(f_mhz) * cable_2m.s21(f_mhz)**2
    cable_2m_c47     = c47.s11(f_mhz) * cable_2m.s21(f_mhz)**2
    cable_2m_c66     = c66.s11(f_mhz) * cable_2m.s21(f_mhz)**2
    cable_0p9m_open  = o.s11(f_mhz) * cable_0p9m.s21(f_mhz)**2
    cable_0p9m_short = s.s11(f_mhz) * cable_0p9m.s21(f_mhz)**2
    cable_0p9m_load  = l.s11(f_mhz) * cable_0p9m.s21(f_mhz)**2
    cable_0p9m_c47   = c47.s11(f_mhz) * cable_0p9m.s21(f_mhz)**2
    cable_0p9m_c66   = c66.s11(f_mhz) * cable_0p9m.s21(f_mhz)**2

    # Compute reflection coefficients for these
    GM_lna  = s2p_lna.s11(f_mhz)
    GM_hot  = s2p_hot.s11(f_mhz)
    GM_cold = s2p_cold.s11(f_mhz)

    # Load antenna S11 -- required for resultant noisewave
    GM_ant = s2p_ant.s11(f_mhz)

    # Generate 1XN matrices for power measured and Gamma
    P_meas  = np.row_stack((#P_2m_open.d(f_mhz),
                            #P_2m_short.d(f_mhz),
                            #P_2m_c47.d(f_mhz),
                            #P_2m_c66.d(f_mhz),
                            #P_cold.d(f_mhz),
                            P_0p9m_open.d(f_mhz),
                            P_0p9m_short.d(f_mhz),
                            P_0p9m_c47.d(f_mhz),
                            #P_0p9m_c66.d(f_mhz),
                            P_cold.d(f_mhz)
                           ))

    GM_meas = np.row_stack((#cable_2m_open,
                            #cable_2m_short,
                            #cable_2m_c47,
                            #cable_2m_c66,
                            #cable_2m_load,
                            cable_0p9m_open,
                            cable_0p9m_short,
                            cable_0p9m_c47,
                            #cable_0p9m_c66,
                            cable_0p9m_load
                           ))

    #############################
    ##   COMPUTE NOISE PARAMS  ##
    #############################

    # Compute noise coefficients C
    C = compute_noise_coeffs(f_mhz,
                   P_hot.d(f_mhz), P_cold.d(f_mhz),
                   T_hot, T_cold,
                   GM_hot, GM_cold, GM_lna,
                   GM_meas, P_meas)

    # Now convert to noise parameters
    NP = C_to_nparams(C)
    R_N, B_opt, G_opt, F_min = NP
    T_min = NF_to_NT(F_min)
    NP[3] = T_min ## Plot temp not F units
    NP_cols = np.column_stack(NP)
    NP_cols = np.column_stack((f_mhz, NP_cols))
    np.savetxt('output_NPs.txt', NP_cols, header='F_MHZ    R_N    B_OPT    G_OPT    T_MIN')

    # And for comparison compute T_rx via Y-factor method
    T_rx_yfm = y_factor(T_hot, T_cold, P_hot.d(f_mhz), P_cold.d(f_mhz))

    # And finally, Compute noisewave for antenna in Kelvin units
    T_nw = compute_T_nw(F_min, R_N, B_opt, G_opt, GM_ant)

    #######
    ### TESTING
    #######

    # Compute FE diode state temperatures
    T_rx_cal = compute_T_nw(F_min, R_N, B_opt, G_opt, (GM_hot+GM_cold)/2)
    #T_rx_cal = T_rx_yfm

    Y_fe_cold = P_fe_cold.d(f_mhz) / P_cold.d(f_mhz)
    T_fe_cold = (Y_fe_cold - 1) * T_rx_cal + Y_fe_cold * T_cold

    Y_fe_hot = P_fe_hot.d(f_mhz) / P_cold.d(f_mhz)
    T_fe_hot = (Y_fe_hot - 1) * T_rx_cal + Y_fe_hot * T_cold


    # 3SS stuff
    P_3ss = (P_cold.d(f_mhz) - P_fe_cold.d(f_mhz)) / (P_fe_hot.d(f_mhz) - P_fe_cold.d(f_mhz))

    T_nw_cal0 = compute_T_nw(F_min, R_N, B_opt, G_opt, GM_cold)
    T_calibrated0 = apply_calibration(P_3ss, T_fe_hot, T_fe_cold,
                      T_rx_cal, T_nw_cal0, GM_cold, GM_lna)


    P_3ss_hot = (P_hot.d(f_mhz) - P_fe_cold.d(f_mhz)) / (P_fe_hot.d(f_mhz) - P_fe_cold.d(f_mhz))
    T_nw_cal_hot = compute_T_nw(F_min, R_N, B_opt, G_opt, GM_hot)
    T_calibrated1 = apply_calibration(P_3ss_hot, T_fe_hot, T_fe_cold,
                      T_rx_cal, T_nw_cal_hot, GM_hot, GM_lna)

    # Difference in Y-factor vs Edward method
    #plt.plot(T_rx_cal / T_rx_yfm)
    #plt.show()

    #plt.figure("3SS test")
    #plt.subplot(211)
    #plt.plot(f_mhz, T_cold)
    #plt.plot(f_mhz, T_calibrated0)
    #plt.subplot(212)
    #plt.plot(f_mhz, T_hot)
    #plt.plot(f_mhz, T_calibrated1)

    #plt.figure("3SS ratios")
    #plt.subplot(211)
    #plt.plot(f_mhz, T_calibrated0 / T_cold)
    #plt.subplot(212)
    #plt.plot(f_mhz, T_calibrated1 / T_hot)


    ###########################
    ##    PLOTTING ROUTINES  ##
    ###########################

    # Plot FE noise diode temps
    plt.figure("Noise diode temps")
    plt.plot(f_mhz, T_fe_cold)
    plt.plot(f_mhz, T_fe_hot)

    # Plot noise parameters
    plt.figure("NPARAM", figsize=(10,8))
    for ii, nparam in enumerate(("R_N", "B_opt", "G_opt", "F_min")):
        plt.subplot(2, 2, ii+1)
        plt.plot(f_mhz, NP[ii])
        plt.title(nparam)
    plt.tight_layout()

    plt.subplot(2,2,1)
    plt.ylim(20, 50)
    plt.subplot(2,2,3)
    plt.ylim(0.010, 0.025)
    plt.subplot(2,2,4)
    plt.ylim(250, 1000)

    for ii in (3,4):
        plt.subplot(2,2,ii)
        plt.xlabel("Frequency [MHz]")
    plt.savefig("NPARAM.png")

    # Plot transducer gain factor term G_S
    plt.figure("G_S")
    G_S = compute_G_S(GM_ant, GM_lna)
    plt.plot(f_mhz, G_S)


    # Plot noise parameters
    plt.figure("NOISEWAVE", figsize=(10,8))
    plt.subplot(2,1,1)
    plt.plot(f_mhz, T_nw,     label='T_rx for antenna')
    plt.plot(f_mhz, T_rx_yfm, label='T_rx by Y-factor')
    plt.legend()
    plt.ylim(300, 1500)
    plt.subplot(2,1,2)
    plt.plot(f_mhz, T_rx_yfm - G_S * T_nw, label="$T_{rx}^{cal} - G_S^{sky} T_{rx}^{sky}$")
    #plt.plot(f_mhz, T_cold, label="$T_{amb}$ ")
    plt.ylim(0, 500)
    plt.legend()

    for ii in (1,2):
        plt.subplot(2,1,ii)
        plt.ylabel("Temperature [K]")
        plt.xlabel("Frequency [MHz]")
    plt.savefig("NOISEWAVE.png")
    plt.show()
