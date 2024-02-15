import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import numpy as np
import matplotlib.pyplot as plt
import math

ydiv = 0.2
wavelength = 632.8e-9
c = 3e8

angleData = csvreader.readData("angle-data.csv")

delta = angleData[0:, 1] / angleData[0:, 2] #equasion 14
ddelta = np.sqrt((0.1 / angleData[0:, 2])**2 + (angleData[0:, 1] / angleData[0:, 2]**2 * 0.1)**2) #using derivatives
nu1 = 2 * np.sqrt(delta) / (1 + delta) #equasion 5
dnu1 = (1 - delta) / np.sqrt(delta) / (1 + delta)**2 * ddelta #using derivatives
nu = (angleData[0:, 4] - angleData[0:, 3]) / (angleData[0:, 4] + angleData[0:, 3]) #equasion 15
dnu = np.sqrt((2 * angleData[0:, 3] / (angleData[0:, 4] + angleData[0:, 3])**2 * 0.1)**2 + (2 * angleData[0:, 4] / (angleData[0:, 4] + angleData[0:, 3])**2 * 0.1)**2) #using derivatives
nu3 = nu / nu1 #equasion 17
dnu3 = np.sqrt((dnu / nu1)**2 + (nu / nu1 * dnu1)**2)

graphs.plot(np.cos(np.radians(angleData[0:, 0])), nu3, np.sin(np.radians(angleData[0:, 0])) * np.radians(2), dnu3, title = "ν3(cosβ)", xlabel = "cosβ", ylabel = "ν3") #laser is likely polarized
graphs.plot(np.cos(np.radians(angleData[0:, 0]))**2, nu3, np.abs(2 * np.cos(np.radians(angleData[0:, 0])) * np.sin(np.radians(angleData[0:, 0])) * np.radians(2)), dnu3, title = "ν3(cos^2β)", xlabel = "cos^2β", ylabel = "ν3")

fig, ax = graphs.basePlot()
plt.title("nu3(beta)")
plt.xlabel("beta, deg")
plt.ylabel("nu3")
ax.errorbar(angleData[0:, 0], nu3, dnu3, 1)
x = np.linspace(angleData[0, 0], angleData[-1, 0], num = 100)
ax.plot(x, np.cos(np.radians(x)), label = "cos(beta)")
ax.plot(x, np.cos(np.radians(x))**2, label = "cos**2(beta)")
plt.legend()
plt.show()

linData = csvreader.readData("lin-data.csv")

delta = linData[0:, 1] / linData[0:, 2] #equasion 14
ddelta = np.sqrt((0.1 / linData[0:, 2])**2 + (linData[0:, 1] / linData[0:, 2]**2 * 0.1)**2) #using derivatives
nu1 = 2 * np.sqrt(delta) / (1 + delta) #equasion 5
dnu1 = (1 - delta) / np.sqrt(delta) / (1 + delta)**2 * ddelta #using derivatives
nu = (linData[0:, 4] - linData[0:, 3]) / (linData[0:, 4] + linData[0:, 3]) #equasion 15
dnu = np.sqrt((2 * linData[0:, 3] / (linData[0:, 4] + linData[0:, 3])**2 * 0.1)**2 + (2 * linData[0:, 4] / (linData[0:, 4] + linData[0:, 3])**2 * 0.1)**2) #using derivatives
nu2 = nu / nu1 #equasion 16
dnu2 = np.sqrt((dnu / nu1)**2 + (nu / nu1 * dnu1)**2)

max1 = 17 #left maximum
dmax1 = 2.5
max2 = 83 #right maximum
dmax2 = 2.5
lhalf = 28 - max1 #l_1/2
dlhalf = 1 + dmax1

L = (max2 - max1) #from graph 1
dL = (dmax1 + dmax2)

print(f"L = ({L:.1f} +- {dL:.1f}) cm")
print(f"Δν_m = ({c / 2 / (L/100) * 1e-6:.0f} +- {c / 2 / (L/100)**2 * dL/100 * 1e-6:.0f}) MHz") #equasion 1

print(f"Δν_full = ({0.6 * c / (lhalf / 100) * 1e-9:.2f} +- {0.6 * c / (lhalf / 100)**2 * (dlhalf / 100) * 1e-9:.2f}) Ghz")
print(f"n = {1 + 0.5 * L / lhalf:.2f} +- {math.sqrt((0.5 * dL / lhalf)**2 + (0.5 * L / lhalf**2 * dlhalf)**2):.2f}")

approx = np.polynomial.polynomial.Polynomial.fit(linData[0:, 0], nu2, 8) #approximate 8 degree polynomial

fig, ax = graphs.basePlot()
plt.title("ν2(x)")
plt.xlabel("x, cm")
plt.ylabel("ν2")
ax.errorbar(linData[0:, 0], nu2, dnu2, 0.1, '.')
x = np.linspace(linData[0, 0], linData[-1, 0], num = 1000)
y = np.linspace(np.min(nu2), np.max(nu2), num = 2)
ax.plot([max1, max1], y, "--r")
ax.plot([max2, max2], y, "--r")
ax.plot(x, approx(x))
plt.show()
