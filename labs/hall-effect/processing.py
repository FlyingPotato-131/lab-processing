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
plt.title('Calibration curve, p3 polynomial')
plt.xlabel('I, A')
plt.ylabel('Φ, Wb')
ax.plot(t, calCurve(t), 'r')
ax.errorbar(calData[0:, 0], calData[0:, 1] / 7.5E-3, 1.5E-4 / 7.5E-3, 0.01, '.')
plt.show()

data = csvreader.readTable("data.csv", 12, titleSize = 0)

fig, ax = plt.subplots()
plt.minorticks_on()
plt.grid(True, "major", "both", color = "#888888")
plt.grid(True, "minor", "both", linestyle = '--')
plt.title('Hall voltage over induction')
plt.xlabel('B')
plt.ylabel('U, V')

slopes = np.empty([0, 2])
for row in data[1:-1]:
	k, b, dk, db = graphs.lsqm(calCurve(data[0, 2:]), row[1] - row[2:], dx = np.full(np.size(data[0, 2:]), dB), dy = np.full(np.size(data[0, 2:]), 1E-6))
	slopes = np.append(slopes, [[k, dk]], axis = 0)
	ax.plot(np.linspace(calCurve(data[0, 2]), calCurve(data[0, -1]), num = 1000), k * np.linspace(calCurve(data[0, 2]), calCurve(data[0, -1]), num = 1000) + b)
	ax.errorbar(calCurve(data[0, 2:]), row[1] - row[2:], np.full(np.size(data[0, 2:]), 1E-6), np.full(np.size(data[0, 2:]), dB), 'b.')
plt.show()
# print(slopes)

k, b, dk, db = graphs.plotLsqm(data[1:-1, 0], slopes[0:, 0], np.full(np.size(data[1: -1, 0]), 1E-5), slopes[0:, 1])

with open('results.txt', 'w') as results:
	results.write('R_h = (' + str(k / 1.5E-3) + '+-' + str(dk / 1.5E-3) + ')m^3 / Cl\n')
	results.write('n = (' + str(1 / k * 1.5E-3 / 1.9E-19 / 1E+18) + '+-' + str(1 / k**2 * 1.5E-6**2 / 1.9E-19 * dk / 1.5E-6 / 1E+18) + ')um^-3\n')
