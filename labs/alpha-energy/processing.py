import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, pi

data = csvreader.readData("data3.csv")

fig, ax = graphs.basePlot()
ax.errorbar(730 - data[:, 0], data[:, 2], 1, 1, fmt = '.')
k, b, dk, db = graphs.lsqm(730 - data[2:10, 0], data[2:10, 2], np.ones(9), np.ones(9))
ax.plot([0, 200], [b, b + 200 * k])
ax.plot([0, 200], [data[1, 2] / 2, data[1, 2] / 2], '--r')

avgP = (data[1, 2] / 2 - b) / k
davgP = sqrt(((0.5 + db) / k)**2 + ((data[1, 2] / 2 - b) / k**2 * dk)**2)
avgx = 9 * avgP / 760
davgx = 9 * davgP / 760
print(f"средний пробег <R> = ({avgx:.3f} +- {davgx:.3f}) см")
rho= 760 * 13546 * 0.01 * 9.8 * 29 / 8.31 / 300 / 1000
print(f"<R> = ({avgx * rho:.3f} +- {davgx * rho:.3f}) г/см^2")
print(f"E = ({(avgx / 0.32)**(2/3):.2f} +- {(davgx / 0.32)**(2/3):.2f}) MeV")

exP = -b / k
dexP = sqrt((db / k)**2 + (db / k**2 * dk)**2)
exx = 9 * exP / 760
dexx = 9 * dexP / 760
print(f"экстаполированный пробег R_ex = ({exx:.3f} +- {dexx:.3f}) см")
print(f"R_ex = ({exx * rho:.2f} +- {dexx * rho:.2f}) г/см^2")
print(f"E = ({(exx / 0.32)**(2/3):.2f} +- {(dexx / 0.32)**(2/3):.2f}) MeV")

plt.title("Зависимость количества частиц от давления")
plt.xlabel("P, torr")
plt.ylabel("N")

plt.show()

fig, ax = graphs.basePlot()
ax.errorbar(730 - data[:, 0], data[:, 1], 1, 1, fmt = '.')
k1, b1, dk1, db1 = graphs.lsqm(730 - data[0:16, 0], data[0:16, 1], np.ones(17), np.ones(17), bflag = False)
k2, b2, dk2, db2 = graphs.lsqm(730 - data[19:, 0], data[19:, 1], np.ones(10), np.ones(10))
ax.plot([0, 600], [b1, b1 + 600 * k1], '--r')
ax.plot([500, 750], [b2 + 500 * k2, b2 + 750 * k2], '--r')

p0 = (b1 - b2) / (k2 - k1)
dp0 = sqrt(((db1 + db2) / (k2 - k1))**2 + ((b1 - b2) / (k2 - k1)**2 * (dk1 + dk2))**2)
x0 = 10 * p0 / 760
dx0 = 10 * dp0 / 760

print(f"средний пробег <R> = ({x0:.2f} +- {dx0:.2f}) см")
print(f"<R> = ({x0 * rho:.1f} +- {dx0 * rho:.1f}) г/см^2")
print(f"E = ({(x0 / 0.32)**(2/3):.1f} +- {(dx0 / 0.32)**(2/3):.1f}) MeV")

exP = 0.5 * p0
dexP = 0.5 * dp0
exx = 10 * exP / 760
dexx = 10 * dexP / 760
print(f"экстаполированный пробег R_ex = ({exx:.2f} +- {dexx:.2f}) см")
print(f"R_ex = ({exx * rho:.1f} +- {dexx * rho:.1f}) г/см^2")
print(f"E = ({(exx / 0.32)**(2/3):.1f} +- {(dexx / 0.32)**(2/3):.1f}) MeV")

plt.show()

# graphs.plot(730 - data[:, 0], data[:, 1])
