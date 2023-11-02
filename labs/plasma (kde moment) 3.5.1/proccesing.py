import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import numpy as np
from numpy.polynomial import Polynomial as poly
import math
import matplotlib.pyplot as plt

pdata = np.array([row for row in csvreader.readData('plasma-vac.csv') if row[0] > 4.9 or row[1] > 3902]) # коэффициенты надо подбирать под задачу если беды с измерениями
# print(pdata)
# graphs.plot(pdata[:, 1] * 1E-6, pdata[:, 0], plotFmt = '.')
fig, ax = plt.subplots()

plt.minorticks_on()
plt.grid(True, "major", "both", color = "#888888")
plt.grid(True, "minor", "both", linestyle = '--')
plt.title('volt-amp of plasma')
plt.xlabel('I, A')
plt.ylabel('V, V')

vac = poly.fit(pdata[:, 1] * 1E-6, pdata[:, 0], 6)
dvac = poly.deriv(vac)
d2vac = poly.deriv(dvac)
Icrit = calculations.newtonMethod(d2vac, 0, 0.0015, 0.01)

I = np.linspace(np.min(pdata[:, 1] * 1E-6), np.max(pdata[:, 1] * 1E-6))
ax.plot(I, vac(I))
# ax.plot(x, d2vac(x))
ax.errorbar(pdata[:, 1] * 1E-6, pdata[:, 0], fmt = '.')
ax.plot(Icrit, vac(Icrit), 'o')
plt.show()

print('наибольшая Rдиф = {} +- {}'.format(dvac(Icrit), 0))
