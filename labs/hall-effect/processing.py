import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import matplotlib.pyplot as plt
import numpy as np
import math
from numpy.polynomial import Polynomial as poly

calData = csvreader.readTable("calibration.csv", 2)

# graphs.plot(calData[0:, 0], calData[0:, 1])
calCurve = poly.fit(calData[0:, 0], calData[0:, 1] / 7.5E-3, 3)
t = np.linspace(calData[0, 0], calData[-1, 0], num = 1000)

dB = math.sqrt(np.sum((calData[0:, 1] / 7.5E-3 - calCurve(calData[0:, 0]))**2 + (1.5E-4 / 7.5E-3)**2)) / np.size(calData[0:, 1])
# print(dB)

fig, ax = plt.subplots()
plt.minorticks_on()
plt.grid(True, "major", "both", color = "#888888")
plt.grid(True, "minor", "both", linestyle = '--')
plt.title('Аппроксимация зависимости B от I')
plt.xlabel('I, А')
plt.ylabel('B, Тл')
ax.plot(t, calCurve(t), 'r')
ax.errorbar(calData[0:, 0], calData[0:, 1] / 7.5E-3, 1.5E-4 / 7.5E-3, 0.01, '.')
plt.show()

data = csvreader.readTable("data.csv", 12, titleSize = 0)

fig, ax = plt.subplots()
plt.minorticks_on()
plt.grid(True, "major", "both", color = "#888888")
# plt.grid(True, "major", "both")
plt.grid(True, "minor", "both", linestyle = '--')
plt.title('График зависимости ЭДС Холла от вектора магнитной индукции')
plt.xlabel('B, Тл')
plt.ylabel('U, мВ')

slopes = np.empty([0, 2])
for row in data[1:-1]:
	k, b, dk, db = graphs.lsqm(calCurve(data[0, 2:]), row[1] - row[2:], dx = np.full(np.size(data[0, 2:]), dB), dy = np.full(np.size(data[0, 2:]), 1E-6))
	slopes = np.append(slopes, [[k, dk]], axis = 0)
	ax.plot(np.linspace(calCurve(data[0, 2]), calCurve(data[0, -1]), num = 1000), 1000 * k * np.linspace(calCurve(data[0, 2]), calCurve(data[0, -1]), num = 1000) + 1000 * b)
	ax.errorbar(calCurve(data[0, 2:]), 1000 * (row[1] - row[2:]), 1000 * np.full(np.size(data[0, 2:]), 1E-6), np.full(np.size(data[0, 2:]), dB), 'b.')
plt.show()
# print(slopes)

k, b, dk, db = graphs.plotLsqm(1000 * data[1:-1, 0], 1000 * slopes[0:, 0], 1000 * np.full(np.size(data[1: -1, 0]), 1E-5), 1000 * slopes[0:, 1], title = "зависимость К от I", xlabel = "I, мA", ylabel = "K * 10^3")

sigma = 1E-3 * 3E-3 / 1.76E-3 / 1.5E-3 / 1.7E-3
dSigma = math.sqrt((1E-5 * 3E-3 / 1.76E-3 / 1.5E-3 / 1.7E-3)**2 + (1E-3 * 3E-3 / 1.76E-3**2 / 1.5E-3 / 1.7E-3 * 1E-6)**2)
n = 1 / k / 1.5E-3 / 1.6E-19
dn = 1 / k**2 / 1.5e-6 / 1.6e-19 * dk * 1.5E-3

with open('results.txt', 'w') as results:
	results.write('R_h = (' + str(k * 1.5E-3) + '+-' + str(dk * 1.5E-3) + ')m^3 / Cl\n')
	results.write('n = (' + str(n) + '+-' + str(dn) + ') m^-3\n')

print('sigma = ' + str(sigma), dSigma)
print('b =', sigma / 1.6E-19 / n, math.sqrt((dSigma / 1.6E-19 / n)**2 + (sigma / 1.6E-19 / n**2 * dn)**2))

