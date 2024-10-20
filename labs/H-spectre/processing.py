import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import matplotlib.pyplot as plt
import numpy as np
import scipy
from math import sqrt, pi

hg = csvreader.readData("hg-data2.csv")
ne = csvreader.readData("ne-data2.csv")
fig, ax = graphs.basePlot()
ax.plot(np.concatenate((hg[:, 0], ne[:, 0])), np.concatenate((hg[:, 1], ne[:, 1] / 10)), '.')

plt.title("calibration")
plt.xlabel("phi, deg")
plt.ylabel("wavelength, nm")

sp = scipy.interpolate.CubicSpline(np.sort(np.concatenate((hg[:, 0], ne[:, 0]))), np.sort(np.concatenate((hg[:, 1], ne[:, 1] / 10))))
x = np.linspace(360, 2472, 1000)
ax.plot(x, sp(x))

print("H_alpha, beta, gamma, delta = ", sp(2427), sp(1432), sp(795), sp(376), "nm")
print("n_1.0, 1.5, gr = ", sp(2320), sp(2218), sp(1605), "nm")
print("R = ", 9 * 4 / (9-4) * 1e9 / sp(2427), 16 * 4 / (16 - 4) * 1e9 / sp(1432), 25 * 4 / (25-4) * 1e9 / sp(795), 36 * 4 / (36 - 4) * 1e9 / sp(376), "m^-1")
hnu2 = (6.6e-34 * 3e8 * 1e9 / 1.6e-19 / sp(2218) - 6.6e-34 * 3e8 * 1e9 / 1.6e-19 / sp(2320)) / 5
print(f"hnu_2 = {hnu2} eV")
print(f"hnu_el = {0.94 - (hnu2 - 0.027)} ev")

plt.show()
# graphs.plot(np.concatenate((hg[:, 0], ne[:, 0])), np.concatenate((hg[:, 1], ne[:, 1] / 10)))
