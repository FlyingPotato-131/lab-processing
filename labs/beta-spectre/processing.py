import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, pi

data = csvreader.readData("data.csv")

data[:, 1] = data[:, 1] - data[-2, 1]
# data[:, 1] = data[:, 1] - 1.1121

graphs.plot(data[:, 0], data[:, 1], 0.01, 0.001, title = "Количество частиц от силы тока в магнитной линзе", xlabel = "I, A", ylabel = "N/t, 1/c")

convpeak = 3.56
pc = data[1:, 0] / convpeak * 875.25 #c * импульс конверсионного пика, кэВ
dpc = np.ones(30) * 0.01 / convpeak * 875.25
E = np.sqrt(511**2 + pc**2) - 511
dE = pc / E * dpc

graphs.plot(E, data[1:, 1], dE, 0.001, title = "Количество частиц от энергии", xlabel = "E, keV", ylabel = "N/t, 1/c")
print(f"разрешение спектрометра R ~ 4.8")
print()


fig, ax = graphs.basePlot()

ax.errorbar(E, np.sqrt(data[1:, 1]) / pc, np.sqrt((0.001 / 2 / np.sqrt(data[1:, 1]) / pc)**2 + (np.sqrt(data[1:, 1]) / pc**2 * dpc)**2), fmt = '.')
k, b, dk, db = graphs.lsqm(E[9:15], np.sqrt(data[9:15, 1]) / pc[9:15], np.sqrt((0.001 / 2 / np.sqrt(data[9:15, 1]) / pc[9:15])**2 + (np.sqrt(data[9:15, 1]) / pc[9:15]**2 * dpc[9:15])**2))

plt.title("График Ферми Кюри")
plt.xlabel("E, keV")
plt.ylabel("N^1/2 / pc")

ax.plot([150, 800], k * np.array([150, 800]) + b)
Ee = -b / k
dEe = sqrt((db / k)**2 + (b / k**2 * dk)**2)
print()
print(f"E_e = ({Ee:.0f} +- {dEe:.0f})кэВ")

plt.show()
