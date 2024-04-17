import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations
import numpy as np
import matplotlib.pyplot as plt
import math

polData = csvreader.readData("pol-data.csv")
graphs.plot(np.cos(np.radians(polData[:, 0]))**2, polData[:, 1], np.ones(np.shape(polData)[0]) * 2 * np.abs(np.cos(np.radians(polData[:, 0])) * np.sin(np.radians(polData[:, 0]))) * np.radians(1), np.ones(np.shape(polData)[0]), title = "intensity over polarization angle", xlabel = "cos^2(α)", ylabel = "receiver voltage V, uV")

interfData = csvreader.readData("interf-data.csv")
print(f"λ = ({2 * (interfData[np.argmax(interfData[:, 1]), 0] - interfData[np.argmin(interfData[:, 1]), 0])} +- {2 * 0.5:.2f}) mm")

fig, ax = graphs.basePlot()
ax.errorbar(interfData[:, 0], interfData[:, 1], 1, 0.25, ".")
k, b, dk, db = graphs.lsqm(interfData[:10, 0], interfData[:10, 1], 0.25 * np.ones(11), np.ones(11))
ax.plot([np.min(interfData[:10, 0]), np.max(interfData[:10, 0])], [k * np.min(interfData[:10, 0]) + b, k * np.max(interfData[:10, 0]) + b])

plt.title("intensity over mirror position")
plt.xlabel("x, mm")
plt.ylabel("receiver voltage V, uV")
plt.show()

maikData = csvreader.readData("maik-data.csv")
k, b, dk, db = graphs.plotlsqm(maikData[:, 0], maikData[:, 1], np.zeros(np.shape(maikData)[0]), np.ones(np.shape(maikData)[0]) * 0.5, title = "position with max intensity over order of maximum", xlabel = "m", ylabel = "x_m, mm")
print(f"λ = ({k} +- {dk}) mm")

I1 = 27
I2 = 39
mirror2Data = csvreader.readData("mirror2-data.csv")
fig, ax = graphs.basePlot()
ax.errorbar(mirror2Data[:, 0] * 0.01, mirror2Data[:, 1], 1, 0.01, '.', label = "experimental data")
x = np.linspace(np.min(mirror2Data[:, 0]) * 0.01, np.max(mirror2Data[:, 0]) * 0.01, num = 1000)
ax.plot(x, I1 + I2 + 2 * math.sqrt(I1 * I2) * np.cos(4 * math.pi / 4.25 * (x + 1.5)), label = "theor data")
plt.legend()
plt.show()

dx = 2
h = 3.2
print(f"n = {dx / h + 1:.3f} +- {0.01 / h:.3f}")
