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
print('ε_0 = {} +- {} Ом * с'.format(e_0, de_0))

# print('σ = {} +- {} (попадает в порядок)'.format(1 / (b * 1E+6)**0.5 * k**0.5 / math.pi / 4.5E-2 / 1.5E-3 / 4 / math.pi / 1E-7, 0))
print('σ_lfe = ({} +- {}) м/Ом'.format(e_0 * kl**0.5 / math.pi / 4.5E-2 / 1.5E-3 / 4 / math.pi / 1E-7, math.sqrt((de_0 * kl**0.5 / math.pi / 4.5E-2 / 1.5E-3 / 4 / math.pi / 1E-7)**2 + (dkl / 2 / math.sqrt(kl) * e_0 / math.pi / 4.5E-2 / 1.5E-3 / 4 / math.pi / 1E-7)**2)))
print('для меди М3 (которая использовалась в работе) σ = 56.2E+6, попадает в порядок')
print()

lfetan = np.array([row for row in lfe if abs(abs(row[3] / row[4] * math.pi - 0.5 * math.pi) - 0.5 * math.pi) > 0.1 and abs(abs(row[3] / row[4] * math.pi - 0.5 * math.pi) - 1.5 * math.pi) > 0.1])
psi0 = lfetan[:, 3] / lfetan[:, 4] * math.pi - 0.5 * math.pi
# print(psi0)
dpsi0 = ((0.1 / lfetan[:, 4] * math.pi)**2 + (0.1 * lfetan[:, 3] / lfetan[:, 4]**2 * math.pi)**2)**0.5

# kt, bt, dkt, dbt = graphs.lsqm(lfetan[:, 0], np.tan(psi0[:]), np.array([1] *  np.size(psi0[:])), dpsi0[:] / np.cos(psi0[:])**2)
# print(kt)
kt, bt, dkt, dbt = graphs.plotlsqm(lfetan[:, 0], np.tan(psi0), np.ones(np.size(lfetan[:, 0])), dpsi0 / np.cos(psi0)**2, title = 'tan(Ψ) of ν', xlabel = 'ν, Hz', ylabel = 'tan(Ψ)')
print('σ_tan = ({} +- {}) м/Ом'.format(kt / math.pi / 4.5E-2 / 1.5E-3 / (4 * math.pi * 1E-7), 0))
print('примерно тот же результат, хз почему такие беды')
print()
# graphs.plot(hfetan[15:, 0], np.tan(psi0[15:]), np.array([1] *  np.size(psi0[15:])), dpsi0[15:] / np.cos(psi0[15:])**2)

hfe = np.array([row for row in data if row[0] >= 400])

psi = hfe[:, 3] / hfe[:, 4] * math.pi - 0.5 * math.pi
dpsi = ((0.1 / hfe[:, 4] * math.pi)**2 + (0.1 * hfe[:, 3] / hfe[:, 4]**2 * math.pi)**2)**0.5
# print(psi)

kh, bh, dkh, dbh = graphs.plotlsqm(hfe[:, 0]**0.5, psi - math.pi / 4, 1 / 2 / hfe[:, 0]**0.5, dpsi, bflag = False, title = 'Ψ of ν^0.5', xlabel = 'ν^0.5, Hz^0.5', ylabel = 'Ψ - π/4, rad')
# print(kh, dkh)
# print(bh)
print('σ_hfe = ({} +- {}) м/Ом'.format(kh**2 / math.pi / 1.5E-3**2 / 4 / math.pi / 1E-7, 2 * kh * dkh / math.pi / 1.5E-3**2 / 4 / math.pi / 1E-7))
print('это сильно ближе, в первом косяков не вижу так что видимо считать из разности фаз точнее')
print()

coildata = csvreader.readData('coil-data.csv')
graphs.plot(coildata[:, 0], coildata[:, 1], 1, 0.1, title = 'L of ν', xlabel = 'ν, Hz', ylabel = 'L, mH')
Lmin = np.min(coildata[:, 1])
Lmax = np.max(coildata[:, 1])

coildata = coildata[(coildata[:, 1] != Lmin)]
graphs.plot(coildata[:, 0]**2, (Lmax - coildata[:, 1]) / (coildata[:, 1] - Lmin), title = 'f(L) of ν^2', xlabel = 'ν^2, Hz^2', ylabel = '(Lmax - L) / (L - Lmin)')
kc, bc, dkc, dbc = graphs.plotLsqm(coildata[:12, 0]**2, (Lmax - coildata[:12, 1]) / (coildata[:12, 1] - Lmin), bflag = False, title = 'f(L) of ν^2, first 12 points, too far from theory to be useful', xlabel = 'ν^2, Hz^2', ylabel = '(Lmax - L) / (L - Lmin)')
# print(kc, dkc)
print('σ_L = ({} +- {}) м/Ом'.format(kc / math.pi / 1.5E-3 / 4.5E-2 / (4 * math.pi * 1E-7), dkc / math.pi / 1.5E-3 / 4.5E-2 / (4 * math.pi * 1E-7)))
print('данный результат бесполезен из-за явного несходства входных данных с теорией')
print()

Hcoeff = data[:, 2] / data[:, 1] / data[:, 0] / e_0
dHcoeff = ((0.001 / data[:, 1] / data[:, 0] / e_0)**2 + (data[:, 2] * 0.001 / data[:, 1]**2 / data[:, 0] / e_0)**2 + (data[:, 2] * 1 / data[:, 1] / data[:, 0]**2 / e_0)**2 + (data[:, 2] * de_0 / data[:, 1] / data[:, 0] / e_0**2)**2)**0.5

h = 1.5E-3
a = 4.5E-2

x = np.linspace(np.min(np.log(data[:, 0])), np.max(np.log(data[:, 0])), num = 1000)
b = np.array([np.sqrt(math.pi * np.exp(x) * e_0 * kl**0.5 / math.pi / a / h), np.sqrt(math.pi * np.exp(x) * kt / math.pi / a / h), np.sqrt(math.pi * np.exp(x) * kh**2 / math.pi / 1.5E-3**2)])
# b = np.sqrt(math.pi * np.exp(x) * 56.2E+6 * 4 * math.pi * 1E-7)
y = 1 / np.sqrt((np.cosh(b*h) * np.cos(b*h) - b*a/2 * np.cosh(b*h) * np.sin(b*h) + b*a/2 * np.sinh(b*h) * np.cos(b*h))**2 + (b*a/2 * np.sinh(b*h) * np.cos(b*h) + np.sinh(b*h) * np.sin(b*h) + b*a/2 * np.cosh(b*h) * np.sin(b*h))**2)
# print(y)
# graphs.plot(x, y)
# graphs.plot(np.log(data[:, 0]), Hcoeff, 1 / data[:, 0], dHcoeff)
fig, ax = plt.subplots()
plt.minorticks_on()
plt.grid(True, "major", "both", color = "#888888")
plt.grid(True, "minor", "both", linestyle = '--')

plt.title('H/H_0 of ln(ν)')
plt.xlabel('ln(ν)')
plt.ylabel('H/H_0')

for i in range(3):
	ax.plot(x, y[i, :], label = ['σ_lfe', 'σ_tan', 'σ_hfe'][i])
ax.errorbar(np.log(data[:, 0]), Hcoeff, 1 / data[:, 0], dHcoeff, fmt = '.', label = 'experimental')
plt.legend()
plt.show()
print('ирония в том, что первые 2 значения проводимости не сходятся с реальным, но показывают верные результаты ослабления поля')
