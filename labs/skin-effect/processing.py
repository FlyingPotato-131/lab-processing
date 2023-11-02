import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import numpy as np
import math
import matplotlib.pyplot as plt

data = csvreader.readData('experiment-data.csv')
# print(data)
lfe = np.array([row for row in data if row[0] <= 200])
e = lfe[:, 2] / lfe[:, 1] / lfe[:, 0]
de = ((0.001 / lfe[:, 1] / lfe[:, 0])**2 + (lfe[:, 2] * 0.001 / lfe[:, 1]**2 / lfe[:, 0])**2 + (lfe[:, 2] * 1 / lfe[:, 1] / lfe[:, 0]**2)**2)**0.5
dy = 2 * de / e**3
dx = 2 * lfe[:, 0] * 1
# print(e)
kl, bl, dkl, dbl = graphs.plotlsqm(lfe[:, 0]**2 * 1E-6, 1 / e**2 * 1E-6, dx * 1E-6, dy * 1E-6, title = '1/ε^2 of ν^2 at lfe', xlabel = 'ν^2, (kHz)^2', ylabel = '1/ε^2, A^2/(V^2*s^2) * 1E+6')
# print(b, db)
# print(kl, dkl)
# print(k)
e_0 = 1 / (bl * 1E+6)**0.5
de_0 = 0.5 * dbl * 1E+6 / (bl * 1E+6)**1.5
# print('ε_0 = {} +- {}'.format(1 / (b * 1E+6)**0.5, 0.5 * db * 1E+6 / (b * 1E+6)**1.5))
print('ε_0 = {} +- {}'.format(e_0, de_0))

# print('σ = {} +- {} (попадает в порядок)'.format(1 / (b * 1E+6)**0.5 * k**0.5 / math.pi / 4.5E-2 / 1.5E-3 / 4 / math.pi / 1E-7, 0))
print('σ_lfe = ({} +- {}) Ом/м'.format(e_0 * kl**0.5 / math.pi / 4.5E-2 / 1.5E-3 / 4 / math.pi / 1E-7, math.sqrt((de_0 * kl**0.5 / math.pi / 4.5E-2 / 1.5E-3 / 4 / math.pi / 1E-7)**2 + (dkl / 2 / math.sqrt(kl) * e_0 / math.pi / 4.5E-2 / 1.5E-3 / 4 / math.pi / 1E-7)**2)))
print('для меди М3 (которая использовалась в работе) σ = 56.2E+6, попадает в порядок')
print()

hfe = np.array([row for row in data if row[0] >= 400])
psi = hfe[:, 3] / hfe[:, 4] * math.pi - 0.5 * math.pi
dpsi = ((0.1 / hfe[:, 4] * math.pi)**2 + (0.1 * hfe[:, 3] / hfe[:, 4]**2 * math.pi)**2)**0.5
# print(psi)

kh, bh, dkh, dbh = graphs.plotlsqm(hfe[:, 0]**0.5, psi - math.pi / 4, 1 / 2 / hfe[:, 0]**0.5, dpsi, bflag = False, title = 'Ψ of ν^0.5', xlabel = 'ν^0.5, Hz^0.5', ylabel = 'Ψ, rad')
# print(kh, dkh)
# print(bh)
print('σ_hfe = ({} +- {}) Ом/м'.format(kh**2 / math.pi / 1.5E-3**2 / 4 / math.pi / 1E-7, 2 * kh * dkh / math.pi / 1.5E-3**2 / 4 / math.pi / 1E-7))
print('это сильно ближе, в первом косяков не вижу так что видимо считать из разности фаз точнее')

Hcoeff = data[:, 2] / data[:, 1] / data[:, 0] / e_0
dHcoeff = ((0.001 / data[:, 1] / data[:, 0] / e_0)**2 + (data[:, 2] * 0.001 / data[:, 1]**2 / data[:, 0] / e_0)**2 + (data[:, 2] * 1 / data[:, 1] / data[:, 0]**2 / e_0)**2 + (data[:, 2] * de_0 / data[:, 1] / data[:, 0] / e_0**2)**2)**0.5

h = 1.5E-3
a = 4.5E-2

x = np.linspace(np.min(np.log(data[:, 0])), np.max(np.log(data[:, 0])), num = 1000)
b = np.sqrt(math.pi * np.exp(x) * kh**2 / math.pi / 1.5E-3**2)
# b = np.sqrt(math.pi * np.exp(x) * 56.2E+6 * 4 * math.pi * 1E-7)
y = 1 / np.sqrt((np.cosh(b*h) * np.cos(b*h) - b*a/2 * np.cosh(b*h) * np.sin(b*h) + b*a/2 * np.sinh(b*h) * np.cos(b*h))**2 + (b*a/2 * np.sinh(b*h) * np.cos(b*h) + np.sinh(b*h) * np.sin(b*h) + b*a/2 * np.cosh(b*h) * np.sin(b*h))**2)
# graphs.plot(x, y)
# graphs.plot(np.log(data[:, 0]), Hcoeff, 1 / data[:, 0], dHcoeff)
fig, ax = plt.subplots()
plt.minorticks_on()
plt.grid(True, "major", "both", color = "#888888")
plt.grid(True, "minor", "both", linestyle = '--')

plt.title('H/H_0 of ln(ν)')
plt.xlabel('ln(ν)')
plt.ylabel('H/H_0')

ax.plot(x, y, label = 'theoretical', color = 'r')
ax.errorbar(np.log(data[:, 0]), Hcoeff, 1 / data[:, 0], dHcoeff, fmt = '.', label = 'experimental')
plt.show()
