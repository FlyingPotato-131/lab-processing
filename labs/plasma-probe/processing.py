import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import math
import numpy as np
import matplotlib.pyplot as plt
import scipy

E0 = 8.85e-12 #universal constant
e = 1.6e-19 #electron charge
K = 1.38e-23 #boltzmann constant
S = math.pi * 0.3e-3 * 7e-3 #probe area
me = 9.109e-31 #electron mass
p = 39.9966 #pressure

single_data = csvreader.readData("single.csv")
double_data = csvreader.readData("double.csv")

fig, ax = graphs.basePlot()
ax.errorbar(single_data[:, 0], single_data[:, 1] * 1e3, 1e-7, 1e-7 * 1e3, fmt = '.') #I(V), one probe
plt.title("I(V) for 1 probe")
plt.xlabel("V, V")
plt.ylabel("I, mA")

ionpoints = 20
k0, delta, dk0, db0 = graphs.lsqm(single_data[:ionpoints, 0], single_data[:ionpoints, 1], 1e-7 * np.ones(ionpoints), 1e-7 * np.ones(ionpoints))
corrected = single_data[:, 1] - (k0 * single_data[:, 0] + delta) #correction for ionic current at Vf
ax.plot([np.min(single_data[:, 0]), np.max(single_data[:, 0])], [k0 * 1e3 * np.min(single_data[:, 0]) + delta * 1e3, k0 * 1e3 * np.max(single_data[:, 0]) + delta * 1e3])
plt.show()

fig, ax = graphs.basePlot()
plt.title("logarithmic for 1 probe")
plt.xlabel("V, V")
plt.ylabel("ln(I)")

rmpoints = 27
l = np.size(corrected)
ax.errorbar(single_data[-rmpoints:, 0], np.log(corrected[-rmpoints:]), 1e-7 / corrected[-rmpoints:], 1e-7, fmt = '.') #log(I) over V
k, b, dk, db = graphs.lsqm(single_data[-rmpoints+5:, 0], np.log(corrected[-rmpoints+5:]), 1e-7 * np.ones(l - rmpoints + 5), 1e-7 / corrected[-rmpoints+5:]) #approximate linear function
ax.plot([np.min(single_data[-rmpoints:, 0]), np.max(single_data[-rmpoints:, 0])], [k * np.min(single_data[-rmpoints:, 0]) + b, k * np.max(single_data[-rmpoints:, 0]) + b]) #plot linear approximation
# ax.plot([np.min(single_data[-rmpoints:, 0]), np.max(single_data[-rmpoints:, 0])], [np.log(-k0 * 217 - delta), np.log(-k0 * 217 - delta)]) #horizontal line
ax.plot([np.min(single_data[-rmpoints:, 0]), np.max(single_data[-rmpoints:, 0])], [np.log(-k0 * np.min(single_data[-rmpoints:, 0]) - delta), np.log(-k0 * np.max(single_data[-rmpoints:, 0]) - delta)]) #horizontal line
ax.plot(single_data[-rmpoints:, 0], np.log(-k0 * single_data[-rmpoints:, 0] - delta)) #horizontal line

Vf = (np.log(-delta) - b) / k #floating potential
dVf = (-1e-7 / delta + db) / k + (np.log(-delta) - b) / k**2 * dk
print(f"Vf = ({Vf:.1f} +- {dVf:.1f}) V")

Te = 11.6e3 / k #electron temperature, 11.6e3 = e / K
dTe = 11.6e3 / k**2 * dk
print(f"T_e = ({Te:.0f} +- {dTe:.0f}) K")

Vp = Vf - K * Te / e * 0.5 * np.log(math.pi / 2 * me / 1.67e-27) #plasma potential
dVp = dVf - K * dTe / e * 0.5 * np.log(math.pi / 2 * me / 1.67e-27)
print(f"Vp = ({Vp:.1f} +- {dVp:.1f}) V")

ne = -(k0 * Vf + delta) * math.exp(0.5) / e / S * math.sqrt(1.67e-27 * 4 / K / Te) #electron concentration
dne = (dk0 * Vf + db0) * math.exp(0.5) / e / S * math.sqrt(1.67e-27 * 4 / K / Te) - (k0 * Vf + delta) * math.exp(0.5) / e / S * math.sqrt(1.67e-27 * 4 / K) / Te**1.5 / 2 * dTe
print(f"ne = ({ne * 1e-12:.0f} +- {dne * 1e-12:.0f}) 1e3/mm^3")

N = p / K / 400 #helium concentration
print(f"N = {N * 1e-12:.0f} 1e3/mm^3")
# dN = p / K / Te**2 * dTe
# print(f"N = ({N * 1e-12:.0f} +- {dN * 1e-12:.0f}) 1e3/mm^3")

a = ne / N #ionization parameter
da = dne / N
print(f"a = ({a * 1e6:.2f} +- {da * 1e6:.2f}) 1e-6")

rD = math.sqrt(E0) * math.sqrt(K * Te / ne) / e #Debai radius
drD = 0.5 * math.sqrt(E0) * math.sqrt(K / ne / Te) / e * dTe + math.sqrt(E0) * math.sqrt(K * Te / ne) / e / ne**1.5 / 2 * dne
print(f"rD = ({rD * 1e6:.2f} +- {drD * 1e6:.2f}) um")
print()

plt.show()

fig, ax = graphs.basePlot()
plt.title("I(V) for 2 probes")
plt.xlabel("V, V")
plt.ylabel("I, mA")

ax.errorbar(double_data[:, 1], double_data[:, 2] * 1e3, 1e-7, 1e-7 * 1e3, fmt = '.') #plot

k1, b1, dk1, db1 = graphs.lsqm(double_data[-30:, 1], double_data[-30:, 2], 1e-7 * np.ones(np.shape(double_data)[0]), 1e-7 * np.ones(np.shape(double_data)[0])) #linears
k2, b2, dk2, db2 = graphs.lsqm(double_data[:30, 1], double_data[:30, 2], 1e-7 * np.ones(np.shape(double_data)[0]), 1e-7 * np.ones(np.shape(double_data)[0]))
k, b, dk, db = graphs.lsqm(double_data[70:90, 1], double_data[70:90, 2], 1e-7 * np.ones(np.shape(double_data)[0]), 1e-7 * np.ones(np.shape(double_data)[0]))

ax.plot([0, double_data[-1, 1]], [b1 * 1e3, k1 * double_data[-1, 1] * 1e3 + b1 * 1e3]) #plot linears
ax.plot([double_data[0, 1], 0], [k2 * double_data[0, 1] * 1e3 + b2 * 1e3, b2 * 1e3])
ax.plot([double_data[70, 1], double_data[90, 1]], [k * double_data[70, 1] * 1e3 + b * 1e3, k * double_data[90, 1] * 1e3 + b * 1e3])

Te = e / K * (-b1 * b2 / (b1 - b2)) / k #electron temperature
dTe = e / K * (b2**2 / (b1 - b2)**2 / k * db1 + b1**2 / (b1 - b2) / k * db2 - b1 * b2 / (b1 - b2) / k**2 * dk)
print(f"T_e = ({Te:.0f} +- {dTe:.0f}) K")

ne = 2 * b1 / e / S * math.sqrt(1.67e-27 * 4 / K / Te) #electron concentration
dne = 2 * db1 / e / S * math.sqrt(1.67e-27 * 4 / K / Te) + 2 * b1 / e / S * math.sqrt(1.67e-27 * 4 / K) / Te**1.5 / 2 * dTe
print(f"ne = ({ne * 1e-12:.0f} +- {dne * 1e-12:.0f}) 1e3/mm^3")

# N = p / K / Te #helium concentration
# dN = p / K / Te**2 * dTe
# print(f"N = ({N * 1e-12:.0f} +- {dN * 1e-12:.0f}) 1e3/mm^3")

a = ne / N #ionization parameter
da = dne / N
print(f"a = ({a * 1e6:.2f} +- {da * 1e6:.2f}) 1e-6")

rD = math.sqrt(E0) * math.sqrt(K * Te / ne) / e #Debai radius
drD = 0.5 * math.sqrt(E0) * math.sqrt(K / ne / Te) / e * dTe + math.sqrt(E0) * math.sqrt(K * Te / ne) / e / ne**1.5 / 2 * dne
print(f"rD = ({rD * 1e6:.2f} +- {drD * 1e6:.2f}) um")
print()

plt.show()