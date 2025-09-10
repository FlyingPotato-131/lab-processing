import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import matplotlib.pyplot as plt
import numpy as np
import math

def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx

cutoff = 400

# read data

spec20 = csvreader.readData("25090320.DAT", split = ' ', datatype = 'U', titlesize = 0) # read data as is (with blank columns)
spec20_x = 0.1 * spec20[cutoff:, 0].astype('f')[::-1] # extract x and y columns and cast them to float, throw away noise, flip array
noise20 = np.average(spec20[:cutoff, 3].astype('f')) # assume white noise
spec20_y = spec20[cutoff:, 3].astype('f')[::-1] - noise20

spec25 = csvreader.readData("25090325.DAT", split = ' ', datatype = 'U', titlesize = 0)
spec25_x = 0.1 * spec25[cutoff:, 0].astype('f')[::-1]
noise25 = np.average(spec25[:cutoff, 3].astype('f'))
spec25_y = spec25[cutoff:, 3].astype('f')[::-1] - noise25

spec30 = csvreader.readData("25090330.DAT", split = ' ', datatype = 'U', titlesize = 0)
spec30_x = 0.1 * spec30[cutoff:, 0].astype('f')[::-1]
noise30 = np.average(spec30[cutoff-10:cutoff, 3].astype('f')) # 30 mA spectre has trailing elements resulting in a shorter noise profile
spec30_y = spec30[cutoff:, 3].astype('f')[::-1] - noise30

# find peaks (damn this will take me twelve years)

J = np.arange(11, 18) # J values

ind20 = np.array([ # automatically determine indices peak values between manually plotted bounds
	find_nearest(spec20_x, 336.2188) + np.argmax(spec20_y[find_nearest(spec20_x, 336.2188) : find_nearest(spec20_x, 336.2505)]),
	find_nearest(spec20_x, 336.1217) + np.argmax(spec20_y[find_nearest(spec20_x, 336.1217) : find_nearest(spec20_x, 336.1532)]),
	find_nearest(spec20_x, 336.0222) + np.argmax(spec20_y[find_nearest(spec20_x, 336.0222) : find_nearest(spec20_x, 336.0561)]),
	find_nearest(spec20_x, 335.9184) + np.argmax(spec20_y[find_nearest(spec20_x, 335.9184) : find_nearest(spec20_x, 335.9432)]),
	find_nearest(spec20_x, 335.7717) + np.argmax(spec20_y[find_nearest(spec20_x, 335.7717) : find_nearest(spec20_x, 335.8011)]),
	find_nearest(spec20_x, 335.6589) + np.argmax(spec20_y[find_nearest(spec20_x, 335.6589) : find_nearest(spec20_x, 335.6816)]),
	find_nearest(spec20_x, 335.5510) + np.argmax(spec20_y[find_nearest(spec20_x, 335.5510) : find_nearest(spec20_x, 335.5777)])
])
wvl20 = spec20_x[ind20] # wavelength values in peaks
int20 = spec20_y[ind20] # intensity values in peaks

ind25 = np.array([
	find_nearest(spec25_x, 336.2052) + np.argmax(spec25_y[find_nearest(spec25_x, 336.2052) : find_nearest(spec25_x, 336.2392)]),
	find_nearest(spec25_x, 336.1080) + np.argmax(spec25_y[find_nearest(spec25_x, 336.1080) : find_nearest(spec25_x, 336.1375)]),
	find_nearest(spec25_x, 336.0086) + np.argmax(spec25_y[find_nearest(spec25_x, 336.0086) : find_nearest(spec25_x, 336.0426)]),
	find_nearest(spec25_x, 335.9049) + np.argmax(spec25_y[find_nearest(spec25_x, 335.9049) : find_nearest(spec25_x, 335.9274)]),
	find_nearest(spec25_x, 335.7537) + np.argmax(spec25_y[find_nearest(spec25_x, 335.7537) : find_nearest(spec25_x, 335.7898)]),
	find_nearest(spec25_x, 335.6455) + np.argmax(spec25_y[find_nearest(spec25_x, 335.6455) : find_nearest(spec25_x, 335.6726)]),
	find_nearest(spec25_x, 335.5307) + np.argmax(spec25_y[find_nearest(spec25_x, 335.5307) : find_nearest(spec25_x, 335.5668)])
])
wvl25 = spec25_x[ind25]
int25 = spec25_y[ind25]

ind30 = np.array([ 
	find_nearest(spec30_x, 336.2188) + np.argmax(spec30_y[find_nearest(spec30_x, 336.2188) : find_nearest(spec30_x, 336.2505)]),
	find_nearest(spec30_x, 336.1217) + np.argmax(spec30_y[find_nearest(spec30_x, 336.1217) : find_nearest(spec30_x, 336.1532)]),
	find_nearest(spec30_x, 336.0222) + np.argmax(spec30_y[find_nearest(spec30_x, 336.0222) : find_nearest(spec30_x, 336.0561)]),
	find_nearest(spec30_x, 335.9184) + np.argmax(spec30_y[find_nearest(spec30_x, 335.9184) : find_nearest(spec30_x, 335.9432)]),
	find_nearest(spec30_x, 335.7717) + np.argmax(spec30_y[find_nearest(spec30_x, 335.7717) : find_nearest(spec30_x, 335.8011)]),
	find_nearest(spec30_x, 335.6589) + np.argmax(spec30_y[find_nearest(spec30_x, 335.6589) : find_nearest(spec30_x, 335.6816)]),
	find_nearest(spec30_x, 335.5510) + np.argmax(spec30_y[find_nearest(spec30_x, 335.5510) : find_nearest(spec30_x, 335.5777)])
])
wvl30 = spec30_x[ind30]
int30 = spec30_y[ind30]

# plot spectres

lw = 0.7
fig, ax = graphs.basePlot()
ax.plot(spec20_x, spec20_y, '-', markersize = 2.5, label = "20 mA", linewidth = lw)
ax.scatter(wvl20, int20, marker = 'x')
ax.plot(spec25_x, spec25_y, '-', markersize = 2.5, label = "25 mA", linewidth = lw)
ax.scatter(wvl25, int25, marker = 'x')
ax.plot(spec30_x, spec30_y, '-', markersize = 2.5, label = "30 mA", linewidth = lw)
ax.scatter(wvl30, int30, marker = 'x')
ax.plot([335.5, 335.5], [0, 5], 'r--', linewidth = 1.5 * lw)
ax.plot([336.0, 336.0], [0, 5], 'r--', linewidth = 1.5 * lw)
ax.plot([336.5, 336.5], [0, 5], 'r--', linewidth = 1.5 * lw)
plt.xlabel("wavelength, nm")
plt.ylabel("intensity, relative units")
plt.title("plasma spectres at different current")
plt.legend()
plt.show()

v = 0 # vibration quantum number
l0inv = 87961.2 + (2047.18 * (v + 0.5) - 28.44 * (v + 0.5)**2) - 58443.2 - (1734.38 * (v + 0.5) - 14.55 * (v + 0.5)**2) # eq 21.11 with E_vib from 21.6
wvl_theor = 1 / (l0inv + 1.824 * J * (J + 1) - 1.638 * J * (J - 1)) * 1e7 #theoretical wavelengths in nm

# error = spec20_x[1] - spec20_x[0] # error equals spectrometer resolution
error = wvl20[0] - wvl25[0] # error accounts for spectrometer play

fig, ax = graphs.basePlot()
ax.errorbar(J, wvl20, error, 0, label = "20 mA", fmt = '.')
ax.errorbar(J, wvl25, error, 0, label = "25 mA", fmt = '.')
ax.errorbar(J, wvl30, error, 0, label = "30 mA", fmt = '.')
ax.plot(J, wvl_theor, label = "theory")
plt.legend()
plt.xlabel("J")
plt.ylabel("wavelength, nm")
plt.title("wavelength in relation to J")
plt.show()

fig, ax = graphs.basePlot()

ax.errorbar(J * (J + 1), np.log(int20 * wvl20**4 / J), 4 * error / wvl20, 0, fmt = '.', label = "20 mA") # 21.13, intensity and wavelength in any units because logarithm
k20, b20, dk20, db20 = graphs.lsqm(J * (J + 1), np.log(int20 * wvl20**4 / J))
ax.plot([11 * 12, 17 * 18], [k20 * 11 * 12 + b20, k20 * 17 * 18 + b20], "tab:blue") # linear approximation plot
Trot20 = -1.824 * 6.626 * 3 / 1.38 * 1e-1 / k20 # 21.14
dTrot20 = 1.824 * 6.626 * 3 / 1.38 * 1e-1 / k20**2 * dk20
print(f"Trot = ({Trot20:.1f} +- {dTrot20:.1f}) K, 20 mA")
print(f"Ttr = ({Trot20 * 1.998 / 1.824:.1f} +- {dTrot20 * 1.998 / 1.824:.1f} K, 20mA)") # 21.15
print()

ax.errorbar(J * (J + 1), np.log(int25 * wvl25**4 / J), 4 * error / wvl25, 0, fmt = '.', label = "25 mA")
k25, b25, dk25, db25 = graphs.lsqm(J * (J + 1), np.log(int25 * wvl25**4 / J))
ax.plot([11 * 12, 17 * 18], [k25 * 11 * 12 + b25, k25 * 17 * 18 + b25], "tab:orange")
Trot25 = -1.824 * 6.626 * 3 / 1.38 * 1e-1 / k25
dTrot25 = 1.824 * 6.626 * 3 / 1.38 * 1e-1 / k25**2 * dk25
print(f"Trot = ({Trot25:.1f} +- {dTrot25:.1f}) K, 25 mA")
print(f"Ttr = ({Trot25 * 1.998 / 1.824:.1f} +- {dTrot25 * 1.998 / 1.824:.1f} K, 25mA)")
print()

ax.errorbar(J * (J + 1), np.log(int30 * wvl30**4 / J), 4 * error / wvl30, 0, fmt = '.', label = "30 mA")
k30, b30, dk30, db30 = graphs.lsqm(J * (J + 1), np.log(int30 * wvl30**4 / J))
ax.plot([11 * 12, 17 * 18], [k30 * 11 * 12 + b30, k30 * 17 * 18 + b30], "tab:green")
Trot30 = -1.824 * 6.626 * 3 / 1.38 * 1e-1 / k30
dTrot30 = 1.824 * 6.626 * 3 / 1.38 * 1e-1 / k30**2 * dk30
print(f"Trot = ({Trot30:.1f} +- {dTrot30:.1f}) K, 30 mA")
print(f"Ttr = ({Trot30 * 1.998 / 1.824:.1f} +- {dTrot30 * 1.998 / 1.824:.1f} K, 30mA)")

plt.xlabel("J(J+1)")
plt.ylabel("ln(I λ^4 / J)")
plt.title("equation 21.13 for determining temperature")

plt.legend()
plt.show()
