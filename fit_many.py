import hickle
import sys
import numpy as np
from scipy.interpolate import interp1d as interp
import scipy.signal
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from lmfit import minimize, Parameters, fit_report


def flag_data(f, d, bp, thr=5):
    """ Flag data. Returns compressed arrays

    flags data where abs(data - bandpass) > threshold

    f: frequency
    d: data (1D)
    bp: bandpass estimate
    thr: Threhold from bandpass above which to flag
    """
    r = d - bp
    d = np.ma.array(d)
    fm = np.ma.array(f)
    d.mask = np.abs(r) > thr
    fm.mask = d.mask
    d.mask[0] = False
    d.mask[-1] = False
    fm.mask[0] = False
    fm.mask[-1] = False
    ff, dd = fm.compressed(), d.compressed()

    return interp(ff, dd)(f)


def fit_poly(x, y, n=5, log=True):		# From 05 script
    """ Fit a polynomial to x, y data 
    
    x (np.array): x-axis of data (e.g. frequency)
    y (np.array): y-axis of data (e.g temperature)
    n (int): number of terms in polynomial (defaults to 5)
    """
    
    x_g = x 
    x = np.ma.array(x, mask=y.mask).compressed()
    y = y.compressed()
    
    if log:
        yl = np.log10(y)
    else:
        yl = y
    fit = np.polyfit(x, yl, n)
    print fit
    p = np.poly1d(fit)
    
    if log:
        return 10**(p(x_g))
    else:
        return p(x_g)

def residual(params, x, model, data):
    mm = model(x, params)
    return (data - mm)

def fit_model_sin_off(x, data):
    params = Parameters()
    params.add('PHI', value=0.3398, vary=True)
    params.add('A_c', value=146., vary=True)
    params.add('PHI0', value=-1.44)
    params.add('B', value=226)
    params.add('M', value=0.2)
    out = minimize(residual, params, args=(x, model_sin_off, data))
    outvals = out.params
    for param, val in out.params.items():
        print "%08s: %2.4f" % (param, val)
    return outvals


def model_sin_off(x, params):
    PHI = params['PHI'].value
    PHI0 = params['PHI0'].value
    A_c = params['A_c'].value
    B = params['B'].value
    M = params['M'].value

    mm = A_c * np.sin(PHI * x + PHI0) + B + M * x
    return mm

def fit_model_damped_sin(x, data):
    params = Parameters()
    params.add("b", value=np.abs(data[0]), vary=True)
    params.add("c", value=1.0, vary=True)
    params.add("d", value=0.0, vary=True)
    params.add("e", value=1.0, vary=True)
    params.add("f", value=0.0, vary=True)
    params.add("g", value=0.0, vary=True)
    out = minimize(residual, params, args=(x, model_damped, data))
    outvals = out.params
    for param, val in out.params.items():
        print "%08s: %2.4f" % (param, val)
    return outvals

def model_damped(x, params):
    b = params["b"].value
    c = params["c"].value
    d = params["d"].value
    e = params["e"].value
    f = params["f"].value
    g = params["g"].value
    mm = (e**(-x))*b*np.sin(c*x+d+f/x)+g
    return mm


lst_min = 11
lst_max = 12
ant = "254A"

accumulated_data = []

for f in sys.argv[1:]:
  data = hickle.load(f)
 
  for key in sorted(data.keys()): 
    if key == "frequencies": frequencies = data[key]
    if key == "lsts": 
      use_lst_indexes = np.arange(data[key].shape[0])[np.logical_and(data[key]>=lst_min, data[key]<=lst_max)]

    if key == ant:
      ant_data = data[key]

  ant_data = np.ma.array(ant_data[use_lst_indexes], mask=ant_data.mask[use_lst_indexes])

  if len(accumulated_data) > 0:
    accumulated_data = np.ma.append(accumulated_data, ant_data, axis=0)
  else: accumulated_data = ant_data

X, Y = np.meshgrid([i for i in range(accumulated_data.shape[1])], [i for i in range(accumulated_data.shape[0])])

print np.ma.min(accumulated_data), np.max(accumulated_data), float(accumulated_data.count())/accumulated_data.shape[0]/accumulated_data.shape[1], "unflagged"
print accumulated_data.shape

ax = plt.figure().gca(projection='3d')
ax.azim = 18
ax.elev = 25
plt.xlabel("Channel")
plt.ylabel("Different times")
ax.plot_surface(X,Y,np.ma.filled(accumulated_data, 0), cmap="rainbow")
plt.show()



# The rest cobbled from leda_cal2/04 script	-------------------------------
D0 = accumulated_data
ff = frequencies


aD = np.ma.mean(D0, axis=0)
rD = aD - fit_poly(ff, aD)

#f2, rD = trim2(rD, ff, 58, 80)	#1167 2084

f2 = ff[1167:2085]
rD = rD[1167:2085]

f2 = np.ma.array(f2, mask=rD.mask).compressed()
rD = rD.compressed()

# Fit a sine wave for flagging
rD_model_params = fit_model_sin_off(f2, rD)
rD_sin_model    = model_sin_off(f2, rD_model_params)
rD = flag_data(f2, rD, rD_sin_model, thr=8)


rD_model_params = fit_model_sin_off(f2, rD)
rD_sin_model    = model_sin_off(f2, rD_model_params)


plt.figure("254A")
plt.subplot(2,1,1)
plt.plot(f2, rD, c='#cc0000', linewidth=0.5)
plt.plot(f2, rD_sin_model, c='#333333')

psin = fit_model_damped_sin(f2, rD - rD_sin_model)
d = model_damped(f2, psin)
res = rD - rD_sin_model
plt.subplot(2,1,2)
plt.plot(f2, rD - rD_sin_model, c='#cc0000')
plt.plot(f2, d, c='#333333', linewidth=0.5)
#plt.savefig("img/04_r254A.png")
plt.show()

filt = res-scipy.signal.medfilt(res, 9)
filt = filt[9:-9]
print np.std(filt)

