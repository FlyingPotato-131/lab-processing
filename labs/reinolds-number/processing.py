import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations
import numpy as np
import matplotlib.pyplot as plt

data = csvreader.readData('data.csv', 0)
# print(data)
Prelation = data[1, 0:] / data[2, 0:]
dPrelation = np.sqrt((1 / data[2, 0:])**2 + (data[1, 0:] / data[2, 0:]**2)**2)
print("ΔP_s / ΔP_core = ", *np.round(Prelation, 3))
print("Δ(ΔP_s / ΔP_core) = ", *np.round(dPrelation, 3))
print()
Urelation = np.sqrt(Prelation)
dUrelation = 0.5 / np.sqrt(Prelation) * dPrelation
print("U / U_core = ", *np.round(Urelation, 3))
print("Δ(U / U_core) = ", *np.round(dUrelation, 3))
print()
Ucore = 1.32 * np.sqrt(data[2, 0:])
dUcore = 1.32 * 0.5 / np.sqrt(data[2, 0:])
print("U_core = ", *np.round(Ucore, 2), " Pa^0.5")
print("Δ(U_core) = ", *np.round(dUcore, 2), " Pa^0.5")
print()
Re = 0.4021e5 * np.sqrt(data[2, 0:])
dRe = 0.4021e5 * 0.5 / np.sqrt(data[2, 0:])
print("Re = ", *np.round(Re / 1000, 1))
print("ΔRe = ", *np.round(dRe / 1000, 1))
print()
print(Re[np.argmax(Urelation)])
print()

fig, ax = graphs.basePlot()
plt.title("U/U_core of Re")
plt.xlabel("Re * 10^3")
plt.ylabel("U/U_core")

ax.plot(Re[np.argmax(Urelation)] / 1000, np.max(Urelation), 'ro')
ax.errorbar(Re / 1000, Urelation, dUrelation, dRe / 1000, '.')
ax.plot(Re / 1000, Urelation, '--')

plt.show()

data110 = csvreader.readData('data-110v.csv', 0)
data220 = csvreader.readData('data-220v.csv', 0)

Url1 = np.sqrt(data110[1, 0:] / data110[1, -1])
dUrl1 = 0.5 / np.sqrt(data110[1, 0:] * data110[1, -1])
Url2 = np.sqrt(data220[1, 0:] / data220[1, -1])
dUrl2 = 0.5 / np.sqrt(data220[1, 0:] * data220[1, -1])

Re1 = 0.4021e5 * np.sqrt(data110[1, 0:])
dRe1 = 0.4021e5 * 0.5 / np.sqrt(data110[1, 0:])
Re2 = 0.4021e5 * np.sqrt(data220[1, 0:])
dRe2 = 0.4021e5 * 0.5 / np.sqrt(data220[1, 0:])

print("U / U_core = ", *np.round(Url1, 3))
print("Re = ", *np.round(Re1 / 1000, 1))
print()
print("U / U_core = ", *np.round(Url2, 3))
print("Re = ", *np.round(Re2 / 1000, 1))
print()

fig, ax = graphs.basePlot()
plt.title("U / Ucore of y")
plt.xlabel("y, mm")
plt.ylabel("U / Ucore")

ax.errorbar(data110[0, 0:], Url1, dUrl1, 0.25, '--xr', label = "110v, laminar")
ax.errorbar(data220[0, 0:], Url2, dUrl2, 0.25, '--.b', label = "220v, turbulent")

plt.legend()
plt.show()
