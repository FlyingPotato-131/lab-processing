import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import numpy as np
from numpy.polynomial import Polynomial as poly
import math
import matplotlib.pyplot as plt
import os
import re
import math
import scipy

def probevac(x, a, b, c):
	return a * np.tanh(b * x) + c * x

pdata = np.array([row for row in csvreader.readData('plasma-vac.csv') if row[0] > 4.9 or row[1] > 3902]) # коэффициенты надо подбирать под задачу если беды с измерениями
# print(pdata)
# graphs.plot(pdata[:, 1] * 1E-6, pdata[:, 0], plotFmt = '.')
fig, ax = plt.subplots()

plt.minorticks_on()
plt.grid(True, "major", "both", color = "#888888")
plt.grid(True, "minor", "both", linestyle = '--')
plt.title('volt-amp of plasma')
plt.ylabel('I, A')
plt.xlabel('U, V')

vac = poly.fit(pdata[:, 1] * 1E-6, pdata[:, 0], 6)
dvac = poly.deriv(vac)
R2 = np.mean(np.sqrt(np.array([(dvac(0.5 * pdata[i, 1] * 1E-6 + 0.5 * pdata[i + 1, 1] * 1E-6) - (pdata[i + 1, 0] - pdata[i, 0]) / (pdata[i + 1, 1] * 1E-6 - pdata[i, 1] * 1E-6))**2 for i in range(np.size(pdata[:, 1]) - 1)])))
# R2 = np.array([(dvac(0.5 * pdata[i, 1] + 0.5 * pdata[i + 1, 1]) - (pdata[i + 1, 0] - pdata[i, 0]) / (pdata[i + 1, 1] - pdata[i, 1]))**2 for i in range(np.size(pdata[:, 1]) - 1)])
print(R2)
d2vac = poly.deriv(dvac)
Icrit = calculations.newtonMethod(d2vac, 0, 0.0015, 0.01)

I = np.linspace(np.min(pdata[:, 1] * 1E-6), np.max(pdata[:, 1] * 1E-6))
ax.plot(vac(I), I)
# ax.plot(x, d2vac(x))
ax.errorbar(pdata[:, 0], pdata[:, 1] * 1E-6, 1E-6, 0.1, fmt = '.')
ax.plot(vac(Icrit), Icrit, 'o')
plt.show()

print('наибольшая Rдиф = ({} +- {}) Ом'.format(dvac(Icrit), R2))
print('участок соответствует участку Ж-З на теоретическом графике')
print()

fig, ax = plt.subplots()

plt.minorticks_on()
plt.grid(True, "major", "both", color = "#888888")
plt.grid(True, "minor", "both", linestyle = '--')
plt.title('volt-amp of probes')
plt.ylabel('I, A')
plt.xlabel('U, V')

# I_m = np.empty([0, 2])
# Rinv = np.empty([0, 2])
T_e = np.empty([0, 2])
ne = np.empty([0, 2])
Ip = np.empty(0)

# print(I_m)
# print(Rinv)

for path in os.listdir('./probe-data/'):
	print(path)
	Ip = np.append(Ip, float(re.search(r'([0-9]*[.])?[0-9]+', path)[0]))

	probe = csvreader.readData('./probe-data/' + path)
	U = np.append(probe[:, 0] - (0.5 * probe[:, 0] - 0.5 * probe[:, 2]), -probe[:, 2] - (0.5 * probe[:, 0] - 0.5 * probe[:, 2]))
	I = np.append(probe[:, 1] - (0.5 * probe[:, 1] - 0.5 * probe[:, 3]), -probe[:, 3] - (0.5 * probe[:, 1] - 0.5 * probe[:, 3])) * 1E-6

	U_0 = np.linspace(np.min(U), np.max(U), num = 1000)
	coeffs, cov = scipy.optimize.curve_fit(probevac, U, I)

	dcoeffs = np.sqrt(np.diag(cov))
	# print(R2)

	ax.plot([np.min(0), np.max(U)], [coeffs[0], coeffs[2] * np.max(U) + coeffs[0]], '--')
	ax.plot(U_0, probevac(U_0, *coeffs))
	ax.plot(U, I, '.')

	print('Iiн = ({} +- {}) A'.format(coeffs[0], dcoeffs[0]))

	Rinv = [coeffs[0] * coeffs[1] + coeffs[2], math.sqrt((dcoeffs[0] * coeffs[1])**2 + (coeffs[0] * dcoeffs[1])**2 + (dcoeffs[2])**2)]
	print('dI/dU = ({} +- {}) Ом^-1'.format(*Rinv))

	#переводим в СГС
	# U = U * 1E+8 / 3E+10
	# I = I * 0.1 * 3E+10
	Iion = coeffs[0] * 0.1 * 3E+10
	dIion = dcoeffs[0] * 0.1 * 3E+10
	Rinv = [r / 1E+9 * 3E+10**2 for r in Rinv]

	T = [0.5 * 4.8E-10 * Iion / Rinv[0] / 1.38E-16, math.sqrt((0.5 * 4.8E-10 * dIion / Rinv[0] / 1.38E-16)**2 + (0.5 * 4.8E-10 * Iion * Rinv[1] / Rinv[0]**2 / 1.38E-16)**2)]
	T_e = np.append(T_e, [[t * 8.62E-5 for t in T]], axis = 0)
	print('Te = ({} +- {}) эВ'.format(*([t * 8.62E-5 for t in T])))
	# print('Te = ({} +- {}) эВ'.format(*([t for t in T])))

	n = [Iion / 0.4 / 4.8E-10 / (math.pi * 0.2E-1 * 5.2E-1) * math.sqrt(20 * 1.672E-24 / 2 / 1.38E-16 / T[0]), math.sqrt((dIion / 0.4 / 4.8E-10 / (math.pi * 0.2E-1 * 5.2E-1) * math.sqrt(20 * 1.672E-24 / 2 / 1.38E-16 / T[0]))**2 + (0.5 * Iion / 0.4 / 4.8E-10 / (math.pi * 0.2E-1 * 5.2E-1) * math.sqrt(20 * 1.672E-24 / 2 / 1.38E-16 / T[0]**2))**2)]
	ne = np.append(ne, [n], axis = 0)
	print('ne = ni = ({} +- {}) cm^-3'.format(*n))

	wp = [math.sqrt(4 * math.pi * n[0] * 4.8E-10**2 / 9.109E-28), 0.5 * math.sqrt(4 * math.pi * 4.8E-10**2 / 9.109E-28 / n[0]) * n[1]]
	print('wp = ({} +- {}) Hz'.format(*wp))

	rDe = [math.sqrt(1.38E-16 * T[0] / 4 / math.pi / n[0] / 4.8E-10**2), math.sqrt((0.5 * math.sqrt(1.38E-16 / 4 / math.pi / n[0] / 4.8E-10**2 / T[0]) * T[1])**2 + (0.5 * math.sqrt(1.38E-16 * T[0] / 4 / math.pi / n[0]**3 / 4.8E-10**2) * n[1])**2)]
	rDi = [math.sqrt(1.38E-16 * 293 / 4 / math.pi / n[0] / 4.8E-10**2), math.sqrt((0.5 * math.sqrt(1.38E-16 / 4 / math.pi / n[0] / 4.8E-10**2 / 293) * 15)**2 + (0.5 * math.sqrt(1.38E-16 * 293 / 4 / math.pi / n[0]**3 / 4.8E-10**2) * n[1])**2)]
	print('rDe = ({} +- {}) cm'.format(*rDe))
	print('rDi = ({} +- {}) cm'.format(*rDi))

	ND = [4/3 * math.pi * n[0] * rDi[0]**3, math.sqrt((4/3 * math.pi * n[1] * rDi[0]**3)**2 + (4 * math.pi * n[0] * rDi[0]**2 * rDi[1])**2)]
	print('ND = {} +- {}'.format(*ND))
	print()

plt.show()
# print(Ip)

# print(T_e)
# print()
# print(ne)
# graphs.plot(Ip, T_e[:, 0], [1E-6] * np.size(Ip), T_e[:, 1], title = 'T_e of Ip', xlabel = 'Ip, A', ylabel = 'T_e, eV')
graphs.plot(Ip, T_e[:, 0], dy = T_e[:, 1], title = 'T_e of Ip', xlabel = 'Ip, mA', ylabel = 'T_e, eV')
graphs.plot(Ip, ne[:, 0], dy = ne[:, 1], title = 'ne of Ip', xlabel = 'Ip, mA', ylabel = 'ne, m^-3') 

# Ip = np.array([])
