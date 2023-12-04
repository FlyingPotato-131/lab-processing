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
h = 0.03
b = 0.005
I = h * b**3 / 12
E = 2e11

data = csvreader.readData('data.csv')

fig, ax = graphs.basePlot()

ky, by, dky, dby = graphs.lsqm(data[0:, 0], data[0:, 1], bflag = False)
kx1, bx1, dkx1, dbx1 = graphs.lsqm(data[0:, 0], data[0:, 2], bflag = False)
kx2, bx2, dkx2, dbx2 = graphs.lsqm(data[0:, 0], data[0:, 3], bflag = False)

ax.plot([np.min(data[0:, 0]), np.max(data[0:, 0])], [ky * np.min(data[0:, 0]), ky * np.max(data[0:, 0])], label = 'Δy approx')
ax.plot([np.min(data[0:, 0]), np.max(data[0:, 0])], [kx1 * np.min(data[0:, 0]), kx1 * np.max(data[0:, 0])], label = 'Δx1 approx')
ax.plot([np.min(data[0:, 0]), np.max(data[0:, 0])], [kx2 * np.min(data[0:, 0]), kx2 * np.max(data[0:, 0])], label = 'Δx2 approx')

ax.plot(data[0:, 0], data[0:, 1], '.', label = 'Δy exp')
ax.plot(data[0:, 0], data[0:, 2], '.', label = 'Δx1 exp')
ax.plot(data[0:, 0], data[0:, 3], '.', label = 'Δx2 exp')

ax.plot(data[0:, 0], data[0:, 0] / 2 * L**3 / 6 / sqrt(2) / E / I * 1e4 * 2, label = 'Δy theor')
ax.plot(data[0:, 0], data[0:, 0] / 2 * L**3 / 6 / sqrt(2) / E / I * 1e4, label = 'Δx theor')

plt.title('offset of P, square frame')
plt.xlabel('P, n')
plt.ylabel('offset, 0.1mm')

plt.legend()
plt.show()

fig, ax = graphs.basePlot()

ky, by, dky, dby = graphs.lsqm(data[0:, 0], data[0:, 5], bflag = False)
kx1, bx1, dkx1, dbx1 = graphs.lsqm(data[0:, 0], data[0:, 6], bflag = False)
kx2, bx2, dkx2, dbx2 = graphs.lsqm(data[0:, 0], data[0:, 7], bflag = False)

ax.plot([np.min(data[0:, 0]), np.max(data[0:, 0])], [ky * np.min(data[0:, 0]), ky * np.max(data[0:, 0])], label = 'Δy approx')
ax.plot([np.min(data[0:, 0]), np.max(data[0:, 0])], [kx1 * np.min(data[0:, 0]), kx1 * np.max(data[0:, 0])], label = 'Δx1 approx')
ax.plot([np.min(data[0:, 0]), np.max(data[0:, 0])], [kx2 * np.min(data[0:, 0]), kx2 * np.max(data[0:, 0])], label = 'Δx2 approx')

ax.plot(data[0:, 0], data[0:, 5], '.', label = 'Δy')
ax.plot(data[0:, 0], data[0:, 6], '.', label = 'Δx1')
ax.plot(data[0:, 0], data[0:, 7], '.', label = 'Δx2')

ax.plot(data[0:, 0], data[0:, 0] / 2 * (pi - 8 / pi) / 4 / E / I * R**3 * 1e4 * 2, label = 'Δy theor')
ax.plot(data[0:, 0], data[0:, 0] / 2 * (8 / pi - 1) / 4 / E / I * R**3 * 1e4, label = 'Δx theor')

plt.title('offset of P, round frame')
plt.xlabel('P, n')
plt.ylabel('offset, 0.1mm')

plt.legend()
plt.show()
