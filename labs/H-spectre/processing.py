import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import matplotlib.pyplot as plt
import numpy as np
import scipy
from math import sqrt, pi

hg = csvreader.readData("hg-data2.csv")
ne = csvreader.readData("ne-data2.csv")
fig, ax = graphs.basePlot()
ax.plot(np.concatenate((hg[:, 0], ne[:, 0])), np.concatenate((hg[:, 1], ne[:, 1] / 10)), '.')

sp = scipy.interpolate.CubicSpline(np.sort(np.concatenate((hg[:, 0], ne[:, 0]))), np.sort(np.concatenate((hg[:, 1], ne[:, 1] / 10))))
x = np.linspace(360, 2472, 1000)
ax.plot(x, sp(x))

print(sp(2427), sp(1432), sp(795), sp(376))
print(sp(2320), sp(2218), sp(1605))

plt.show()
# graphs.plot(np.concatenate((hg[:, 0], ne[:, 0])), np.concatenate((hg[:, 1], ne[:, 1] / 10)))
