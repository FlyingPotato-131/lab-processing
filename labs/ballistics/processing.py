import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations
import matplotlib.pyplot as plt
import numpy as np
import math

def sign(x):
	return 1 if x > 0 else -1 if x < 0 else 0


L = 12e-3
Lb = 0.5
D = 6.85e-3
m = 0.6e-3
pdiv = 0.4 #atm per mark

gamma = 1.4
c = 340
a2 = gamma / (2 * c**2)
S = math.pi * D**2 / 4

data = csvreader.readData("data.csv")

# pressure = np.linspace(np.min(data[0:, 0]) * pdiv, np.max(data[0:, 0]) * pdiv, num = 100)
pressure = np.linspace(1, np.max(data[0:, 0]) * pdiv, num = 100)
N = 10000
velocity = np.empty(100)
dx = Lb / N
for n in range(100):
	tmpv2 = 0
	for i in range(N):
		du2 = S / m / 2 * (pressure[n] * 1e5 * abs(1 - (gamma - 1) / 2 * tmpv2 / c**2)**(gamma / (gamma - 1)) - (1 + (gamma - 1) / 2 * tmpv2 / c**2)**(gamma / (gamma - 1)) * 1e5) * dx
		tmpv2 = tmpv2 + du2
	velocity[n] = math.sqrt(abs(tmpv2))

print(f"u_inf0 = {1 / math.sqrt(a2):.0f} m/s (theoretical velocity after infinite acceleration with 0 ext pressure)")
print(f"u_inf = {1 / math.sqrt(a2) * math.sqrt((20 - 1) / (20 + 1)):.0f} m/s (theoretical velocity after infinite acceleration with 1 atm ext pressure, 20 atm int pressure)")

fig, ax = graphs.basePlot()
plt.title("speed over pressure")
plt.xlabel("p, atm")
plt.ylabel("v, m/s")

ax.errorbar(data[0:, 0] * pdiv, L / data[0:, 1] * 1e6, 0.5 * pdiv, np.sqrt((0.1e-3 / data[0:, 1])**2 + (L / data[0:, 1]**2 * 0.1e6)**2), '.', label = "experimental data")
ax.plot(pressure, 1 / math.sqrt(a2) * np.sqrt(1 - np.exp(-2 * S * pressure * 1e5 / m * a2 * Lb)), label = "theor, 0 ext pressure")
ax.plot(pressure, 1 / math.sqrt(a2) * np.sqrt((pressure - 1) / (pressure + 1)) * (1 - np.exp(-2 * S * 1e5 / m * a2 * (pressure + 1) * Lb)), label = "theor, ext atmospheric pressure")
ax.plot(pressure, velocity, label = "numerical integration")

plt.legend()
plt.show()
