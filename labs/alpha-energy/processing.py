import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, pi

data1 = csvreader.readData("data3.csv")

graphs.plot(730 - data1[:, 0], data1[:, 1])
# graphs.plot(730 - data1[:, 0], data1[:, 2])
# data1.sort(axis = 0)
# data1 = data1[::-1]


fig, ax = graphs.basePlot()
ax.errorbar(730 - data1[:, 0], data1[:, 2], fmt = '.')
k, b, dk, db = graphs.lsqm(730 - data1[2:10, 0], data1[2:10, 2])
ax.plot([0, 200], [b, b + 200 * k])

print(k, b)

plt.show()
