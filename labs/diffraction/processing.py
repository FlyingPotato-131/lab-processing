import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations
import numpy as np
import math
import matplotlib.pyplot as plt

data1 = csvreader.readData("data1.csv")
k, b, dk, db = graphs.plotlsqm(data1[:, 3], np.sin(np.radians(data1[:, 0] + data1[:, 1] / 60)), np.zeros(np.shape(data1)[0]), np.cos(np.radians(data1[:, 0] + data1[:, 1] / 60)) * 1 / 60, title = "diffraction angle over wavelength", ylabel = "sin(φ_m)", xlabel = "λ, nm", bflag = False)

print(f"d = ({1 / k * 1e-3:.4f} +- {1 / k**2 * dk * 1e-3:.4f}) μm")

data2 = csvreader.readData("data2.csv")
dwvl = 20 #angstrem
wvl = 5780
Dright = -(data2[:, 2] * 60 + data2[:, 3] - data2[:, 5] * 60 - data2[:, 6]) / dwvl
Dleft = (data2[:, 7] * 3600 + data2[:, 8] * 60 + data2[:, 9] - data2[:, 10] * 3600 - data2[:, 11] * 60 - data2[:, 12]) / dwvl
Dtheor = data2[:, 0] / np.sqrt(1 / k**2 * 100 - data2[:, 0]**2 * wvl**2) * 180 / math.pi * 3600
x = np.linspace(-np.max(data2[:, 0]), np.max(data2[:, 0]), num = 1000)
DtheorGraph = abs(x / np.sqrt(1 / k**2 * 100 - x**2 * wvl**2) * 180 / math.pi * 3600)
print(f"Dright = {Dright}")
print(f"Dleft = {Dleft}")
np.set_printoptions(precision=2)
print(f"Dtheor = {Dtheor}")

fig, ax = graphs.basePlot()
plt.title("D over m")
plt.xlabel("m")
plt.ylabel("D, sec/angstrem")
ax.scatter(data2[:, 0], Dright)
ax.scatter(-data2[:, 0], Dleft)
ax.plot(x, DtheorGraph)
plt.show()

data3 = csvreader.readData("data3.csv")
phi1 = data3[:, 1] + data3[:, 2] / 60 + data3[:, 3] / 3600
phi2 = data3[:, 4] + data3[:, 5] / 60 + data3[:, 6] / 3600
wvl1 = 1 / k * np.sin(np.radians(phi1)) / data3[:, 0]
wvl2 = 1 / k * np.sin(np.radians(phi2)) / data3[:, 0]
R = wvl1 / np.abs(wvl1 - wvl2)
print(f"R = {R}")