import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations
import math
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial as poly

data = csvreader.readTable("data.csv", 3)
misc = csvreader.readMisc("misc-data.csv", 2, hasTitle = True)

xaxis = np.array([data[i, 0] + (misc[1] * 10**-3) * data[i, 1] + 273 for i in range(np.shape(data)[0])])
yaxis = np.array([10**-12 * 1 / ((value * 10**-6)**2 - (misc[0] * 10**-6)**2) for value in data[0:, 2]])

dx = np.array([math.sqrt((0.01)**2 + (misc[1] * 10**-3 * 1)**2) for i in range(np.shape(data)[0])])
dy = np.array([10**-12 * 2 * value * 0.01 * 10**-12 / ((value * 10**-6)**2 - (misc[0] * 10**-6)**2)**2 for value in data[0:, 2]])

fig, ax = plt.subplots()
plt.minorticks_on()
plt.grid(True, "major", "both", color = "#888888")
plt.grid(True, "minor", "both", linestyle = '--')
plt.title("Зависимость 1/χ от Τ")
plt.xlabel("T, K")
plt.ylabel("1/Χ, МГц^2")

k, b, dk, db = graphs.lsqm(xaxis[-4:-1], yaxis[-4:-1], dx, dy)
# Dk = math.sqrt(dk**2 + (dy[-1] / dx[-1])**2 + (yaxis[-1] * dx[-1] / xaxis[-1]**2)**2)	
# Db = math.sqrt(db**2 + dy[-1]**2 + (k * dx[-1])**2)
# print(dk, db)
theta_p = -b / k
dTheta_p = math.sqrt((db / k)**2 + (b * dk / k**2)**2)
ax.plot([-b / k, xaxis[-1]], [0, k * xaxis[-1] + b], 'r')

approx = poly.fit(xaxis[0:-3], yaxis[0:-3], 4)
theta_k = calculations.newtonMethod(approx, 0, 280, 0.001)
index = calculations.closestValue(xaxis, theta_k)
R2 = np.average(np.array([abs(approx(xaxis[i]) - yaxis[i]) + 0.5 * dy[index] + 0.5 * dy[index + 1] for i in range(np.shape(data)[0])]))
dTheta_k = R2 / approx.deriv()(theta_k)
t = np.linspace(xaxis[0], xaxis[-4], num = 1000)
ax.plot(t, approx(t), 'r')

ax.errorbar(xaxis, yaxis, yerr = dy, xerr = dx, fmt='.')

# graphs.plotPoly(7, xaxis, yaxis, dx = dx, dy = dy)

with open('results.txt', 'w') as results:
	results.write('Θ_p = (' + str(theta_p) + ' +- ' + str(dTheta_p) + ')' + '\n')
	results.write('Θ_k = (' + str(theta_k) + ' +- ' + str(dTheta_k) + ')' + '\n')

plt.savefig('graph.svg')
plt.show()
