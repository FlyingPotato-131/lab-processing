import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations
import math
import matplotlib.pyplot as plt
import numpy as np

wavelength = 546e-9 #possibly incorrect

coordsReal = [0, 0.8]
coords = [-328, 865]
delta = (coordsReal[1] - coordsReal[0]) / (coords[1] - coords[0]) # цена деления, мм/дел
ddelta = math.sqrt((2 * 0.01 / (coords[1] - coords[0]))**2 + (2 * (coordsReal[1] - coordsReal[0]) / (coords[1] - coords[0])**2 * 1)**2)

data = csvreader.readData("data.csv", 0)

kdark, bdark, dkdark, dbdark = graphs.lsqm(data[2:, 0], data[2:, 1] * delta, np.zeros(data.shape[0]), np.sqrt((delta)**2 + (data[0:, 1] * ddelta)**2), bflag = 1) #МНК без учета первых 2 точек ибо там искажения от механических сдвигов
klight, blight, dklight, dblight = graphs.lsqm(data[2:, 0], data[2:, 2] * delta, np.zeros(data.shape[0]), np.sqrt((delta)**2 + (data[0:, 2] * ddelta)**2), bflag = 1) #МНК без учета первых 2 точек ибо там искажения от механических сдвигов

Rdark = kdark / wavelength
dRdark = dkdark / wavelength

Rlight = klight / wavelength
dRlight = dklight / wavelength

print(f"R1 = ({1e-3 * Rdark:.1f} +- {1e-3 * dRdark:.1f}) m")
print(f"R2 = ({1e-3 * Rlight:.1f} +- {1e-3 * dRlight:.1f}) m")

fig, ax = graphs.basePlot()

plt.title("squared radii over ring number")
plt.xlabel("m")
plt.ylabel("r^2, mm^2")

ax.errorbar(data[0:, 0], data[0:, 1] * delta, np.sqrt((delta)**2 + (data[0:, 1] * ddelta)**2), 0, '.', label = "dark rings") #графики из п.3
ax.plot([data[0, 0], data[-1, 0]], [kdark * data[0, 0] + bdark, kdark * data[-1, 0] + bdark], label = "dark approx")

ax.errorbar(data[0:, 0], data[0:, 2] * delta, np.sqrt((delta)**2 + (data[0:, 2] * ddelta)**2), 0, '.', label = "light rings")
ax.plot([data[0, 0], data[-1, 0]], [klight * data[0, 0] + blight, klight * data[-1, 0] + blight], label = "light approx")

plt.legend()
plt.show()
