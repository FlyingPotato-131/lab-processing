import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import math
import matplotlib.pyplot as plt
import numpy as np

nring = csvreader.readData("nring.csv")

graphs.plotlsqm(nring[:, 1] / 100, nring[:, 0], title = "n rings at frequency", xlabel = "freq, Hz", ylabel = "n rings")

nring2 = csvreader.readData("nring2.csv")

graphs.plotlsqm(nring2[:, 1] / 100, nring2[:, 0], title = "n rings at frequency, outer ~1Hz", xlabel = "freq, Hz", ylabel = "n rings")

Re1 = [2055, 1935, 1760, 1575, 1390, 1185, 975, 785, 632, 833, 1206, 1645, 2111, 2590, 3076, 3565, 4057] #граница устойчивости из теории
Re2 = [-4000, -3500, -3000, -2500, -2000, -1500, -1000, -500, 0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000]

stab = csvreader.readData("stab.csv") #граница устойчивости
dstab0 = 1 * 3e-2**2 / 2e-5 / 100 * 2 * math.pi + stab[:, 0] * 3e-2**2 / 2e-5**2 / 100 * 2 * math.pi * 0.05 * 2e-5
dstab1 = 1 * 3.5e-2**2 / 2e-5 / 100 * 2 * math.pi + stab[:, 0] * 3.5e-2**2 / 2e-5**2 / 100 * 2 * math.pi * 0.05 * 2e-5
stab[:, 0] = stab[:, 0] * 3e-2**2 / 2e-5 / 100 * 2 * math.pi #пересчет в Re
stab[:, 1] = stab[:, 1] * 3.5e-2**2 / 2e-5 / 100 * 2 * math.pi

negext = -204 * 3.5e-2**2 / 2e-5 / 100 * 2 * math.pi #противонаправленное вращение
dnegext = 1 * 3.5e-2**2 / 2e-5 / 100 * 2 * math.pi + 204 * 3.5e-2**2 / 2e-5**2 / 100 * 2 * math.pi * 0.05 * 2e-5
negint = 760 * 3e-2**2 / 2e-5 / 100 * 2 * math.pi
dnegint = 1 * 3e-2**2 / 2e-5 / 100 * 2 * math.pi + 760 * 3e-2**2 / 2e-5**2 / 100 * 2 * math.pi * 0.05 * 2e-5

fig, ax = graphs.basePlot()

plt.title("stability edge")
plt.xlabel("Re ext")
plt.ylabel("Re int")

ax.plot(Re2, Re1, label = "theory")
ax.errorbar(stab[:, 1], stab[:, 0], dstab0, dstab1, ".", label = "same direction")
ax.errorbar(negext, negint, dnegint, dnegext, ".", label = "reversed")
plt.legend()

plt.show()