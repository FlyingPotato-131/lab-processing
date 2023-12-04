import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import math
import numpy as np
import matplotlib.pyplot as plt

TofCdata = csvreader.readData('TofC-0.1Hn-0Ohm.csv')

L0 = 0.1
Texp = TofCdata[0:, 1] / TofCdata[0:, 2]
dTexp = 0.01e-3 / TofCdata[0:, 2]
Ttheor = 2 * math.pi * np.sqrt(L0 * TofCdata[0:, 0] * 1e-6)
dTtheor = math.pi * np.sqrt(L0 / TofCdata[0:, 0] / 1e-6) * 0.0001e-6

fig, ax = graphs.basePlot();
plt.title('Texp of Ttheor')
plt.xlabel('Ttheor, μs')
plt.ylabel('Texp, μs')
ax.errorbar(Ttheor * 1e6, Texp * 1e6, dTexp * 1e6, dTtheor * 1e6, '.')
plt.show()

thetaofRdata = csvreader.readData('thetaofR-0.1Hn-6nF.csv')

C1 = 6e-9
Theta = 1 / thetaofRdata[0:, 3] * np.log(thetaofRdata[0:, 1] / thetaofRdata[0:, 2])
dTheta = np.sqrt((1 / thetaofRdata[0:, 3] / thetaofRdata[0:, 1] * 1e-3)**2 + (1 / thetaofRdata[0:, 3] / thetaofRdata[0:, 2] * 1e-3)**2)
RL1 = 53 #берем сопротивление, соответствующее ближайшей частоте
R1 = thetaofRdata[0:, 0] + RL1
dR1 = 2

k2, b2, dk2, dbd = graphs.lsqm(np.sort(1 / R1**2)[0:-3], np.sort(1 / Theta**2)[0:-3], np.sort(1 / R1**3 * dR1)[0:-3], np.sort(1 / Theta**3 * dTheta)[0:-3])
print(f'Rcr_exp = {2 * math.pi * math.sqrt(k2):.0f} +- {math.pi * dk2 / math.sqrt(k2):.0f} Ohm')
print(f'Rcr_theor = {2 * math.sqrt(L0 / C1):.0f} Ohm\n')
print(f'Qmax = {math.pi / np.min(Theta):.3f} +- {math.pi / np.min(Theta)**2 * np.min(dTheta):.3f}') #наибольшая добротность при наименьшем коэффициенте затухания
print(f'Qmin = {math.pi / np.max(Theta):.3f} +- {math.pi / np.max(Theta)**2 * np.max(dTheta):.3f}\n') #наоборот

fig, ax = graphs.basePlot()
plt.title('1/Θ^2 of 1/R^2')
plt.xlabel('1/R^2, kOhm^-2')
plt.ylabel('1/Θ^2')
ax.errorbar(1 / R1**2 * 1e6, 1 / Theta**2, 2 / Theta**3 * dTheta, 2 / R1**3 * dR1 * 1e6, '.')
ax.plot([np.min(1 / R1**2 * 1e6), np.max(1 / R1**2 * 1e6)], [k2 * np.min(1 / R1**2), k2 * np.max(1 / R1**2)])
plt.show()
# print(Theta)

thetaofRdata2d = csvreader.readData('thetaofR2D-0.1Hn-6nF.csv')
Theta2d = 1 / thetaofRdata2d[0:, 3] * np.log(thetaofRdata2d[0:, 1] / (thetaofRdata2d[0:, 1] - thetaofRdata2d[0:, 2]))
dTheta2d = np.sqrt((thetaofRdata2d[0:, 2] / (thetaofRdata2d[0:, 1] - thetaofRdata2d[0:, 2]) / thetaofRdata2d[0:, 1] * 0.1)**2 + (0.1 / (thetaofRdata2d[0:, 1] - thetaofRdata2d[0:, 2]))**2 ) / thetaofRdata2d[0:, 3]
print(f'Qmax2d = {math.pi / np.min(Theta2d):.2f} +- {math.pi / np.min(Theta2d)**2 * np.min(dTheta2d):.2f}') #наибольшая добротность при наименьшем коэффициенте затухания
print(f'Qmin2d = {math.pi / np.max(Theta2d):.2f} +- {math.pi / np.max(Theta2d)**2 * np.min(dTheta2d):.2f}\n') #наоборот

print(f'Qmax_theor = {math.sqrt(L0 / C1) / 400:.2f}')
print(f'Qmin_theor = {math.sqrt(L0 / C1) / 2000:.2f}\n')

R3 = 400
R4 = 800

apfc400 = csvreader.readData('apfc3-400Ohm.csv')
apfc800 = csvreader.readData('apfc4-800Ohm.csv')

U0nu0400 = apfc400[np.argmax(apfc400, axis = 0)[1]]
U0nu0800 = apfc800[np.argmax(apfc800, axis = 0)[1]]

# print(apfc400[0:np.shape(apfc400)[0] // 2, 1])
# print(apfc400[0:np.argmax(apfc400, axis = 0)[1], 1])
# print(apfc400[np.argmax(apfc400, axis = 0)[1]:, 1])
# print(apfc400)

dw400 = [calculations.closestValue(apfc400[0:np.argmax(apfc400, axis = 0)[1], 1], U0nu0400[1] / math.sqrt(2)), calculations.closestValue(apfc400[np.argmax(apfc400, axis = 0)[1]:, 1][::-1], U0nu0400[1] / math.sqrt(2))]
# print(dw400)
kinvl400 = (apfc400[dw400[0] + 1, 0] - apfc400[dw400[0], 0]) / (apfc400[dw400[0] + 1, 1] - apfc400[dw400[0], 1]) * U0nu0400[1] / U0nu0400[0]
binvl400 = apfc400[dw400[0], 0] / U0nu0400[0] - kinvl400 * apfc400[dw400[0], 1] / U0nu0400[1]
kinvr400 = (apfc400[-2 - dw400[1], 0] - apfc400[-1 - dw400[1], 0]) / (apfc400[-2 - dw400[1], 1] - apfc400[-1 - dw400[1], 1]) * U0nu0400[1] / U0nu0400[0]
binvr400 = apfc400[-1 - dw400[1], 0] / U0nu0400[0] - kinvr400 * apfc400[-1 - dw400[1], 1] / U0nu0400[1]
dkinv400 = (apfc400[dw400[0] + 1, 0] - apfc400[dw400[0], 0]) / (apfc400[dw400[0] + 1, 1] - apfc400[dw400[0], 1])**2 * 0.02 * U0nu0400[1] / U0nu0400[0]
dbinv400 = dkinv400 * apfc400[dw400[0], 1] / U0nu0400[1]

dw800 = [calculations.closestValue(apfc800[0:np.argmax(apfc800, axis = 0)[1], 1], U0nu0800[1] / math.sqrt(2)), calculations.closestValue(apfc800[np.argmax(apfc800, axis = 0)[1]:, 1][::-1], U0nu0800[1] / math.sqrt(2))]
# print(dw800)
kinvl800 = (apfc800[dw800[0] + 1, 0] - apfc800[dw800[0], 0]) / (apfc800[dw800[0] + 1, 1] - apfc800[dw800[0], 1]) * U0nu0800[1] / U0nu0800[0]
binvl800 = apfc800[dw800[0], 0] / U0nu0800[0] - kinvl800 * apfc800[dw800[0], 1] / U0nu0800[1]
kinvr800 = (apfc800[-2 - dw800[1], 0] - apfc800[-1 - dw800[1], 0]) / (apfc800[-2 - dw800[1], 1] - apfc800[-1 - dw800[1], 1]) * U0nu0800[1] / U0nu0800[0]
binvr800 = apfc800[-1 - dw800[1], 0] / U0nu0800[0] - kinvr800 * apfc800[-1 - dw800[1], 1] / U0nu0800[1]
dkinv800 = (apfc800[dw800[0] + 1, 0] - apfc800[dw800[0], 0]) / (apfc800[dw800[0] + 1, 1] - apfc800[dw800[0], 1])**2 * 0.02 * U0nu0800[1] / U0nu0800[0]
dbinv800 = dkinv800 * apfc800[dw800[0], 1] / U0nu0800[1]

print(f'Q400afc = {apfc400[np.argmax(apfc400, axis = 0)[1], 0] / (kinvr400 / math.sqrt(2) + binvr400 - kinvl400 / math.sqrt(2) - binvl400) / U0nu0400[0]:.2f} +- {apfc400[np.argmax(apfc400, axis = 0)[1], 0] / (kinvr400 / math.sqrt(2) + binvr400 - kinvl400 / math.sqrt(2) - binvl400)**2 * (dkinv400 / math.sqrt(2) + dbinv400) / U0nu0400[0]:.2f}')
print(f'Q800afc = {apfc800[np.argmax(apfc800, axis = 0)[1], 0] / (kinvr800 / math.sqrt(2) + binvr800 - kinvl800 / math.sqrt(2) - binvl800) / U0nu0800[0]:.2f} +- {apfc800[np.argmax(apfc800, axis = 0)[1], 0] / (kinvr800 / math.sqrt(2) + binvr800 - kinvl800 / math.sqrt(2) - binvl800)**2 * (dkinv800 / math.sqrt(2) + dbinv800) / U0nu0800[0]:.2f}\n')

fig, ax = graphs.basePlot()
plt.title('AFC of oscillator')
plt.xlabel('ν/ν0 where ν0 is the resonance frequency')
plt.ylabel('U/U0 where U0 is the voltage at resonance')
ax.plot(apfc400[0:, 0] / U0nu0400[0], apfc400[0:, 1] / U0nu0400[1], '.-', label = '400 Ohm')
ax.plot(apfc800[0:, 0] / U0nu0800[0], apfc800[0:, 1] / U0nu0800[1], '.-', label = '800 Ohm')
# ax.plot([apfc400[dw400[0], 0] / U0nu0400[0], apfc400[1 - dw400[1], 0] / U0nu0400[0]], [apfc400[dw400[0], 1] / U0nu0400[1], apfc400[1 - dw400[1], 1] / U0nu0400[1]])
ax.plot([kinvl400 / math.sqrt(2) + binvl400, kinvr400 / math.sqrt(2) + binvr400], [1 / math.sqrt(2), 1 / math.sqrt(2)], 'o')
ax.plot([kinvl800 / math.sqrt(2) + binvl800, kinvr800 / math.sqrt(2) + binvr800], [1 / math.sqrt(2), 1 / math.sqrt(2)], 'o')
plt.legend()
plt.show()

invapfc400 = np.empty([0, 3])
for row in apfc400:
	if 2 * math.pi * row[2] * row[0] > 0.5 * math.pi:
		invapfc400 = np.append(invapfc400, [row], axis = 0)

apfc400[0:, 2] = -apfc400[0:, 2]
invapfc400[0:, 2] = -invapfc400[0:, 2]

# print(-math.pi - 2 * math.pi * invapfc400[0:, 2] * invapfc400[0:, 0])
dw400 = [calculations.closestValue(-math.pi - 2 * math.pi * invapfc400[0:, 2] * invapfc400[0:, 0], -0.25 * math.pi), calculations.closestValue(2 * math.pi * apfc400[0:, 2][::-1] * apfc400[0:, 0][::-1], -0.25 * math.pi)]

kinvl400 = 2 * math.pi * (invapfc400[dw400[0] + 1, 0] - invapfc400[dw400[0], 0]) / (-2 * math.pi * invapfc400[dw400[0] + 1, 2] * invapfc400[dw400[0] + 1, 0] + 2 * math.pi * invapfc400[dw400[0], 2] * invapfc400[dw400[0], 0])
binvl400 = 2 * math.pi * invapfc400[dw400[0], 0] - kinvl400 * (-math.pi - 2 * math.pi * invapfc400[dw400[0], 2] * invapfc400[dw400[0], 0])
kinvr400 = (apfc400[-2 - dw400[1], 0] - apfc400[-1 - dw400[1], 0]) / (apfc400[-2 - dw400[1], 2] * apfc400[-2 - dw400[1], 0] - apfc400[-1 - dw400[1], 2] * apfc400[-1 - dw400[1], 0])
binvr400 = 2 * math.pi * apfc400[-1 - dw400[1], 0] - 2 * math.pi * kinvr400 * apfc400[-1 - dw400[1], 2] * apfc400[-1 - dw400[1], 0]
dkinv400 = (invapfc400[dw400[0] + 1, 0] - invapfc400[dw400[0], 0]) / (invapfc400[dw400[0] + 1, 2] * invapfc400[dw400[0] + 1, 0] - invapfc400[dw400[0], 2] * invapfc400[dw400[0], 0])**2 * 1e-7 * invapfc400[dw400[0], 0]
dbinv400 = dkinv400 * invapfc400[dw400[0], 2] * invapfc400[dw400[0], 0]

invapfc800 = np.empty([0, 3])
for row in apfc800:
	if 2 * math.pi * row[2] * row[0] > 0.5 * math.pi:
		invapfc800 = np.append(invapfc800, [row], axis = 0)

apfc800[0:, 2] = -apfc800[0:, 2]
invapfc800[0:, 2] = -invapfc800[0:, 2]

dw800 = [calculations.closestValue(-math.pi - 2 * math.pi * invapfc800[0:, 2] * invapfc800[0:, 0], -0.25 * math.pi), calculations.closestValue(2 * math.pi * apfc800[0:, 2][::-1] * apfc800[0:, 0][::-1], -0.25 * math.pi)]

kinvl800 = 2 * math.pi * (invapfc800[dw800[0] + 1, 0] - invapfc800[dw800[0], 0]) / (-2 * math.pi * invapfc800[dw800[0] + 1, 2] * invapfc800[dw800[0] + 1, 0] + 2 * math.pi * invapfc800[dw800[0], 2] * invapfc800[dw800[0], 0])
binvl800 = 2 * math.pi * invapfc800[dw800[0], 0] - kinvl800 * (-math.pi - 2 * math.pi * invapfc800[dw800[0], 2] * invapfc800[dw800[0], 0])
kinvr800 = (apfc800[-2 - dw800[1], 0] - apfc800[-1 - dw800[1], 0]) / (apfc800[-2 - dw800[1], 2] * apfc800[-2 - dw800[1], 0] - apfc800[-1 - dw800[1], 2] * apfc800[-1 - dw800[1], 0])
binvr800 = 2 * math.pi * apfc800[-1 - dw800[1], 0] - 2 * math.pi * kinvr800 * apfc800[-1 - dw800[1], 2] * apfc800[-1 - dw800[1], 0]
dkinv800 = (invapfc800[dw800[0] + 1, 0] - invapfc800[dw800[0], 0]) / (invapfc800[dw800[0] + 1, 2] * invapfc800[dw800[0] + 1, 0] - invapfc800[dw800[0], 2] * invapfc800[dw800[0], 0])**2 * 1e-7 * invapfc800[dw800[0], 0]
dbinv800 = dkinv800 * invapfc800[dw800[0], 2]

print(f'Q400pfc = {2 * math.pi * invapfc400[-1, 0] / (-kinvr400 * 0.25 * math.pi + binvr400 + kinvl400 * 0.25 * math.pi - binvl400):.2f} +- {2 * math.pi * invapfc400[-1, 0] / (-kinvr400 * 0.25 * math.pi + binvr400 + kinvl400 * 0.25 * math.pi - binvl400)**2 * (dkinv400 * 0.25 * math.pi + dbinv400):.2f}')
print(f'Q800pfc = {2 * math.pi * invapfc800[-1, 0] / (-kinvr800 * 0.25 * math.pi + binvr800 + kinvl800 * 0.25 * math.pi - binvl800):.2f} +- {2 * math.pi * invapfc800[-1, 0] / (-kinvr800 * 0.25 * math.pi + binvr800 + kinvl800 * 0.25 * math.pi - binvl800)**2 * (dkinv800 * 0.25 * math.pi + dbinv800):.2f}\n')

fig, ax = graphs.basePlot()
plt.title('PFC of oscillator')
plt.xlabel('ω, Hz')
plt.ylabel('Δφ')
ax.plot([2 * math.pi * np.min(apfc400[0:, 0]), 2 * math.pi * np.max(apfc400[0:, 0])], [-0.5 * math.pi, -0.5 * math.pi], 'r--')

ax.plot(2 * math.pi * apfc400[0:, 0], 2 * math.pi * apfc400[0:, 2] * apfc400[0:, 0], '.-', label = '400 Ohm')
ax.plot(2 * math.pi * apfc800[0:, 0], 2 * math.pi * apfc800[0:, 2] * apfc800[0:, 0], '.-', label = '800 Ohm')

ax.plot(2 * math.pi * invapfc400[0:, 0], -1.0 * math.pi - 2 * math.pi * invapfc400[0:, 2] * invapfc400[0:, 0], '.--', label = '400 Ohm mirror')
ax.plot(2 * math.pi * invapfc800[0:, 0], -1.0 * math.pi - 2 * math.pi * invapfc800[0:, 2] * invapfc800[0:, 0], '.--', label = '800 Ohm mirror')

ax.plot([-kinvl400 * 0.25 * math.pi + binvl400, -kinvr400 * 0.25 * math.pi + binvr400], [-0.25 * math.pi, -0.25 * math.pi], 'o')
ax.plot([-kinvl800 * 0.25 * math.pi + binvl800, -kinvr800 * 0.25 * math.pi + binvr800], [-0.25 * math.pi, -0.25 * math.pi], 'o')
plt.legend()
plt.show()

U0400up = 8.28
U1400up = 7.24e-1
U8400up = 7.96
theta400up = 1 / 7 * math.log((U0400up - U1400up) / (U0400up - U8400up))
dtheta400up = 1 / 7 * (U0400up - U8400up) / (U0400up - U1400up) * 3 * 0.1 / (U0400up - U8400up)
print(f'Q400up = {math.pi / theta400up:.3f} +- {math.pi / theta400up**2 * dtheta400up:.3f}')

U1400lo = 5.6
U7400lo = 6.4e-1
theta400lo = 1 / 6 * math.log(U1400lo / U7400lo)
dtheta400lo = 1 / 6 * U7400lo / U1400lo * 0.1 / U7400lo
print(f'Q400lo = {math.pi / theta400lo:.3f} +- {math.pi / theta400lo**2 * dtheta400lo:.3f}')

U0800up = 8.28
U1800up = 7.24e-1
U8800up = 7.96
theta800up = 1 / 5 * math.log((U0800up - U1800up) / (U0800up - U8800up))
dtheta800up = 1 / 5 * (U0800up - U8800up) / (U0800up - U1800up) * 3 * 0.1 / (U0800up - U8800up)
print(f'Q800up = {math.pi / theta800up:.3f} +- {math.pi / theta800up**2 * dtheta800up:.3f}')

U1800lo = 4.4
U7800lo = 0.3
theta800lo = 1 / 4 * math.log(U1800lo / U7800lo)
dtheta800lo = 1 / 4 * U7800lo / U1800lo * 0.1 / U7800lo
print(f'Q800lo = {math.pi / theta800lo:.3f} +- {math.pi / theta800lo**2 * dtheta800lo:.3f}')
