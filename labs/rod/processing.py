import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import matplotlib.pyplot as plt

data = csvreader.readData('data.csv')

fig, ax = graphs.basePlot()

k, b, dk, db = graphs.lsqm(data[0:, 0], data[0:, 1])

ktheor = 0.03
# надо заменить на значение из теории

ax.plot(data[0:, 0], data[0:, 1], '.')
ax.plot(data[0:, 0], k * data[0:, 0] + b)
ax.plot(data[0:, 0], ktheor * data[0:, 0] + b)

plt.show()
