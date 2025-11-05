import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import matplotlib.pyplot as plt
import numpy as np

D = 60 #mm
dx = 150 #mm
x0 = 70 #mm
dr = 10 #mm (assumption)
densw = 1000 #water density
densa = 1.2 #air density (approximate)
g = 10

data = csvreader.readData("data.csv", titlesize = 0)

np.set_printoptions(linewidth = 400, precision = 0) #make debugging easier

#change units in data matrix
data[:, 0] = data[:, 0] - x0 #set correct zero position for x
data[:, 1:] = np.sqrt(2 * densw / densa * g * (data[0, 1:] - data[:, 1:]) / 1000) #change height to velocity in m/s
# print(data)

#plot velocity profiles
r = np.arange(-150, 151, dr)
fig, ax = graphs.basePlot()
for h in range(data.shape[0] - 1):
	ax.plot(r + 50, data[h + 1, 1:], label = f"x = {data[h + 1, 0]:.0f} mm") #r + 50 to center profiles at zero
plt.title("velocity profiles")
plt.xlabel("r, mm")
plt.ylabel("U, m/s")
plt.legend()
plt.show()

#r = -50 <=> col 11 is flow center

#plot dimensionless velocity over distance
fig, ax = graphs.basePlot()
ax.plot(data[1:, 0] / D, data[1:, 11] / data[1, 11])
plt.title("dimensionless axial velocity")
plt.xlabel("x / D")
plt.ylabel("U / U[x = 0]")
plt.show()

#plot flow width
fig, ax = graphs.basePlot()
widths = np.empty(9)
for h in range(data.shape[0] - 1):
	Uhalf = np.max(data[h + 1, 1:]) / 2
	width = 10 * data[h + 1, 1:][data[h + 1, 1:] > Uhalf].size
	widths[h] = width
ax.plot(data[1:, 0], widths)
plt.title("flow width along x")
plt.xlabel("x, mm")
plt.ylabel("r, mm")
plt.show()

#plot dimensionless velocity profiles
r = np.arange(-150, 151, dr)
fig, ax = graphs.basePlot()
for h in range(data.shape[0] - 1):
	ax.plot((r + 50) / D, data[h + 1, 1:] / data[h + 1, 11], label = f"x = {data[h + 1, 0]:.0f} mm")
plt.title("velocity profiles")
plt.xlabel("r / D")
plt.ylabel("U / U[r = 0]")
plt.legend()
plt.show()
