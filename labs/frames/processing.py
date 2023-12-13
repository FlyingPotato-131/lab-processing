import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, pi

L = 0.3 / sqrt(2)
R = 0.395 / 2
h = 0.0302
b = 0.005
I = h * b**3 / 12
E = 2e11

data = csvreader.readData('data.csv')

fig, ax = graphs.basePlot()

ky, by, dky, dby = graphs.lsqm(data[0:, 0], data[0:, 1], bflag = True)
kx1, bx1, dkx1, dbx1 = graphs.lsqm(data[0:, 0], data[0:, 2], bflag = True)
kx2, bx2, dkx2, dbx2 = graphs.lsqm(data[0:, 0], data[0:, 3], bflag = True)

ax.plot([np.min(data[0:, 0]), np.max(data[0:, 0])], [ky * np.min(data[0:, 0]) + by, ky * np.max(data[0:, 0]) + by], label = 'Δy approx')

ax.plot(data[0:, 0], data[0:, 1], '.', label = 'Δy exp')

ax.plot(data[0:, 0], data[0:, 0] * L**3 / 6 / sqrt(2) / E / I * 1e4 + by, label = 'Δy theor')

plt.title('offset of P, square frame')
plt.xlabel('P, n')
plt.ylabel('offset, 0.1mm')

print(f'ktheor = {L**3 / 6 / sqrt(2) / E / I * 1e4}')
print(f'ky = {ky} +- {dky}\n')

plt.legend()
plt.show()

fig, ax = graphs.basePlot()

ax.plot([np.min(data[0:, 0]), np.max(data[0:, 0])], [kx1 * np.min(data[0:, 0]) + bx1, kx1 * np.max(data[0:, 0]) + bx1], label = 'Δx1 approx')
ax.plot([np.min(data[0:, 0]), np.max(data[0:, 0])], [kx2 * np.min(data[0:, 0]) + bx2, kx2 * np.max(data[0:, 0]) + bx2], label = 'Δx2 approx')

ax.plot(data[0:, 0], data[0:, 2], '.', label = 'Δx1 exp')
ax.plot(data[0:, 0], data[0:, 3], '.', label = 'Δx2 exp')

# ax.plot(data[0:, 0], data[0:, 0] / 2 * L**3 / 6 / sqrt(2) / E / I * 1e4 + bx1, label = 'Δx theor')
ax.plot(data[0:, 0], data[0:, 0] * L**3 / 12 / sqrt(2) / E / I * 1e4 + bx1, label = 'Δx theor')

plt.title('offset of P, square frame')
plt.xlabel('P, n')
plt.ylabel('offset, 0.1mm')

print(f'ktheor = {L**3 / 12 / sqrt(2) / E / I * 1e4}')
print(f'kx = {kx1} +- {dkx1}\n')

plt.legend()
plt.show()

h = 0.0295
b = 0.0055
I = h * b**3 / 12

fig, ax = graphs.basePlot()

ky, by, dky, dby = graphs.lsqm(data[0:, 0], data[0:, 5], bflag = True)
kx1, bx1, dkx1, dbx1 = graphs.lsqm(data[0:, 0], data[0:, 6], bflag = True)
kx2, bx2, dkx2, dbx2 = graphs.lsqm(data[0:, 0], data[0:, 7], bflag = True)

ax.plot([np.min(data[0:, 0]), np.max(data[0:, 0])], [ky * np.min(data[0:, 0]) + by, ky * np.max(data[0:, 0]) + by], label = 'Δy approx')

ax.plot(data[0:, 0], data[0:, 5], '.', label = 'Δy')

# ax.plot(data[0:, 0], data[0:, 0] / 2 * (pi - 8 / pi) / 32 / E / I * R**3 * 1e4 * 2 + by, label = 'Δy theor')
# ax.plot(data[0:, 0], data[0:, 0] / 2 * 2 * R**3 / E / I * (pi / 4 - 2 / pi) * 1e4, label = 'Δy theor')
ax.plot(data[0:, 0], data[0:, 0] * R**3 / E / I * (pi / 4 - 2 / pi) * 1e4, label = 'Δy theor')
# ax.plot(data[0:, 0], data[0:, 0] * )
#пока здесь подгнаны 16 и 32, надо уточнить рассчеты

plt.title('offset of P, round frame')
plt.xlabel('P, n')
plt.ylabel('offset, 0.1mm')

print(f'ktheor = {R**3 / E / I * (pi / 4 - 2 / pi) * 1e4}')
print(f'ky = {ky} +- {dky}')

plt.legend()
plt.show()

fig, ax = graphs.basePlot()

ax.plot([np.min(data[0:, 0]), np.max(data[0:, 0])], [kx1 * np.min(data[0:, 0]) + bx1, kx1 * np.max(data[0:, 0]) + bx1], label = 'Δx1 approx')
ax.plot([np.min(data[0:, 0]), np.max(data[0:, 0])], [kx2 * np.min(data[0:, 0]) + bx2, kx2 * np.max(data[0:, 0]) + bx2], label = 'Δx2 approx')

ax.plot(data[0:, 0], data[0:, 6], '.', label = 'Δx1')
ax.plot(data[0:, 0], data[0:, 7], '.', label = 'Δx2')

# ax.plot(data[0:, 0], data[0:, 0] / 2 * (8 / pi - 1) / 16 / E / I * R**3 * 1e4, label = 'Δx theor')
ax.plot(data[0:, 0], data[0:, 0] * R**3 / 2 / E / I * (2 / pi) * 1e4, label = 'Δx theor')

plt.title('offset of P, round frame')
plt.xlabel('P, n')
plt.ylabel('offset, 0.1mm')

print(f'ktheor = {R**3 / 2 / E / I * (2 / pi) * 1e4}')
print(f'k = {kx1} +- {dkx1}')

plt.legend()
plt.show()
