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
graphs.plot(data1[:, 0], data1[:, 1])

data2 = csvreader.readData("2.62v.csv")
graphs.plot(data2[:, 0], data2[:, 1])
data2[:, 1] = data2[:, 1] + 0.05

data3 = csvreader.readData("2.94v.csv")
graphs.plot(data3[:, 0], data3[:, 1])
data3[:, 1] = data3[:, 1] + 0.05

fig, ax = graphs.basePlot()
ax.errorbar(data1[:, 0], -np.log(data1[:, 1] / np.max(data1[:, 1])), fmt = '.')
ax.errorbar(data2[:, 0], -np.log(data2[:, 1] / np.max(data2[:, 1])), fmt = '.')
ax.errorbar(data3[:, 0], -np.log(data3[:, 1] / np.max(data3[:, 1])), fmt = '.')
plt.show()