import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, pi

data1 = csvreader.readData("2.29v.csv")
data1[:, 1] = data1[:, 1] + 0.05
data1[:, 0] = data1[:, 0] - 1

E2 = data1[np.argmin(data1[np.argmax(data1[:, 1]):, 1]) + np.argmax(data1[:, 1]), 0]
E1 = data1[np.argmax(data1[:, 1]), 0]
print(f"E1 = ({E1:.2f} +- 0.01) ev")
print(f"E2 = ({E2:.2f} +- 0.01) ev")
print(f"l= ({6.6e-34 * sqrt(5) / sqrt(32 * 9.109e-31 * (E2 - E1) * 1.6e-19) * 1e9:.3f} +- {6.6e-34 * sqrt(5) / sqrt(32 * 9.109e-31) / 2 / sqrt((E2 - E1) * 1.6e-19)**3 * 1 * 1.6e-19 * 1e9:.3f}) nm")
print(f"U = ({4/5 * E2 - 9/5 * E1:.3f} +- {4/5 * 0.01 + 9/5 * 0.01:.3f}) ev")
print()

graphs.plot(data1[:, 0], data1[:, 1], np.ones(np.size(data1[:, 0])) * 0.01, np.ones(np.size(data1[:, 0])) * 0.01, title = "зависимость анодного тока от напряжения катод-сетка при V_канал = 2.29V", xlabel = "U_к-с, V", ylabel = "I_а * 100kOhm, mV")

data2 = csvreader.readData("2.62v.csv")

data2[:, 1] = data2[:, 1] + 0.05
data2[:, 0] = data2[:, 0] - 1

E2 = data2[np.argmin(data2[np.argmax(data2[:, 1]):, 1]) + np.argmax(data2[:, 1]), 0]
E1 = data2[np.argmax(data2[:, 1]), 0]
print(f"E1 = ({E1:.2f} +- 0.01) ev")
print(f"E2 = ({E2:.2f} +- 0.01) ev")
print(f"l= ({6.6e-34 * sqrt(5) / sqrt(32 * 9.109e-31 * (E2 - E1) * 1.6e-19) * 1e9:.3f} +- {6.6e-34 * sqrt(5) / sqrt(32 * 9.109e-31) / 2 / sqrt((E2 - E1) * 1.6e-19)**3 * 1 * 1.6e-19 * 1e9:.3f}) nm")
print(f"U = ({4/5 * E2 - 9/5 * E1:.3f} +- {4/5 * 0.01 + 9/5 * 0.01:.3f}) ev")
print()

# graphs.plot(data2[:, 0], data2[:, 1])
graphs.plot(data2[:, 0], data2[:, 1], np.ones(np.size(data2[:, 0])) * 0.01, np.ones(np.size(data2[:, 0])) * 0.01, title = "зависимость анодного тока от напряжения катод-сетка при V_канал = 2.62", xlabel = "U_к-с, V", ylabel = "I_а * 100kOhm, mV")

data3 = csvreader.readData("2.94v.csv")

data3[:, 1] = data3[:, 1] + 0.05
data3[:, 0] = data3[:, 0] - 1

E2 = data3[np.argmin(data3[np.argmax(data3[:, 1]):, 1]) + np.argmax(data3[:, 1]), 0]
E1 = data3[np.argmax(data3[:, 1]), 0]
print(f"E1 = ({E1:.2f} +- 0.01) ev")
print(f"E2 = ({E2:.2f} +- 0.01) ev")
print(f"l= ({6.6e-34 * sqrt(5) / sqrt(32 * 9.109e-31 * (E2 - E1) * 1.6e-19) * 1e9:.3f} +- {6.6e-34 * sqrt(5) / sqrt(32 * 9.109e-31) / 2 / sqrt((E2 - E1) * 1.6e-19)**3 * 1 * 1.6e-19 * 1e9:.3f}) nm")
print(f"U = ({4/5 * E2 - 9/5 * E1:.3f} +- {4/5 * 0.01 + 9/5 * 0.01:.3f}) ev")
print()

# graphs.plot(data3[:, 0], data3[:, 1])
graphs.plot(data3[:, 0], data3[:, 1], np.ones(np.size(data3[:, 0])) * 0.01, np.ones(np.size(data3[:, 0])) * 0.01, title = "зависимость анодного тока от напряжения катод-сетка при V_канал = 2.94", xlabel = "U_к-с, V", ylabel = "I_а * 100kOhm, mV")

fig, ax = graphs.basePlot()
ax.errorbar(data1[3:, 0], -np.log(data1[3:, 1] / np.max(data1[3:, 1])), np.max(data1[3:, 1]) / data1[3:, 1] * 0.01, 0.01, fmt = '.-', label = "2.29v")
ax.errorbar(data2[2:, 0], -np.log(data2[2:, 1] / np.max(data2[2:, 1])), np.max(data2[2:, 1]) / data2[2:, 1] * 0.01, 0.01, fmt = '.-', label = "2.62v")
ax.errorbar(data3[3:, 0], -np.log(data3[3:, 1] / np.max(data3[3:, 1])), np.max(data3[3:, 1]) / data3[3:, 1] * 0.01, 0.01, fmt = '.-', label = "2.94v")

plt.title("логарифмическая зависимость I_анод от напряжения катод-сетка")
plt.xlabel("U_катод-сетка, V")
plt.ylabel("-ln(I_анод / I_max)")
plt.legend()

plt.show()

fig, ax = graphs.basePlot()
ax.errorbar(data1[:, 0], data1[:, 1], np.ones(np.size(data1[:, 0])) * 0.01, np.ones(np.size(data1[:, 0])) * 0.01, fmt = '.-', label = "2.29v")
ax.errorbar(data2[:, 0], data2[:, 1], np.ones(np.size(data2[:, 0])) * 0.01, np.ones(np.size(data2[:, 0])) * 0.01, fmt = '.-', label = "2.62v")
ax.errorbar(data3[:, 0], data3[:, 1], np.ones(np.size(data3[:, 0])) * 0.01, np.ones(np.size(data3[:, 0])) * 0.01, fmt = '.-', label = "2.94v")
plt.title("зависимость анодного тока от напряжения катод-сетка")
plt.xlabel("U_к-с, V")
plt.ylabel("I_а * 100kOhm, mV")
plt.legend()

plt.show()