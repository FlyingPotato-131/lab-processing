import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

from math import pi
import math
import numpy as np

h = 6.6e-34
mu_B = 927.4e-26
N = 49
D = 14.3e-3
w = 50 * 2 * pi
Bcoeff = 1 / N / (pi * D**2 / 4) / w
dBcoeff = 1 / N / (pi * D**2 / 4)**2 / w * pi * D * 0.1e-3 / 2

calib = csvreader.readData("calib.csv")
k, x, dk, dx = graphs.plotlsqm(calib[:, 0], 1e4 * (0.5e-3 * calib[:, 1] + 0.5e-3 * calib[:, 2]) / N / (pi * D**2 / 4) / w, 0.1 * np.ones(np.shape(calib)[0]), 1e4 * np.sqrt((0.1e-3 / N / (pi * D**2 / 4) / w)**2 + (1e-3 * calib[:, 1] / N / (pi * D**2 / 4)**2 / w * pi * D * 0.1e-3 / 2)**2), title = "calibration graph", xlabel = "main coil voltage, mV", ylabel = "base magnetic field, Gs", bflag = False)
print("хз нахер нам этот график нужен, мы без него поле измеряли")

B0 = 11.6e-3 * Bcoeff * 1e4
dB0 = math.sqrt((0.1e-3 * Bcoeff)**2 + (11.6e-3 * dBcoeff)**2) * 1e4
print(f"B0 = ({B0:.2f} +- {dB0:.2f}) Gs")

nu0 = 124.31e6
print(f"ν_0 = (124.31 +- 0.01) MHz")

g = h * nu0 / mu_B / B0 * 1e4
dg = math.sqrt((h * 0.01e6 / mu_B / B0 * 1e4)**2 + (h * nu0 / mu_B / B0**2 * dB0 * 1e4)**2)
print(f"g = {g:.3f} +- {dg:.3f}")

Bmod = 7.21e-3 * Bcoeff * 1e4
dBmod = math.sqrt((0.1e-3 * Bcoeff)**2 + (7.21e-3 * dBcoeff)**2) * 1e4
print(f"Bmod = ({Bmod:.2f} +- {dBmod:.2f}) Gs")
B = 2 * Bmod / 8 * 4/5
dB = 2 * dBmod / 8 * 4/5
print(f"B = ({B:.2f} +- {dB:.2f}) Gs")
print(f"ν = ({2 * mu_B * B / h / 1e9:.1f} +- {2 * mu_B * dB / h / 1e9:.1f}) GHz")
