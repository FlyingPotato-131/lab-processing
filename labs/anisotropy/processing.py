import sys

sys.path.append("../../")
import csvreader
import graphs
import numpy as np
import matplotlib.pyplot as plt
import math

data = csvreader.readTable('data.csv', 4, titleSize = 1)

sigma = data[0:, 2] * 10 / data[0:, 1] * 1000 / data[0:, 3] * 1000

x = np.linspace(0, math.pi / 2, num = 1000)

fig, ax = plt.subplots()

plt.minorticks_on()
plt.grid(True, "major", "both", color = "#888888")
plt.grid(True, "minor", "both", linestyle = '--')

plt.ylabel('нормированный предел прочности')
plt.xlabel('угол, рад')

ax.plot(data[0:, 0] / 180 * math.pi, sigma / sigma[0], '.')

ax.plot(x, [0.58 / math.sqrt(0.58**2 * math.cos(x[i])**4 + math.sin(x[i])**4 + 1.12 * math.sin(x[i])**2 * math.cos(x[i])**2) for i in range(len(x))])

#graphs.plot(data[0:, 0], sigma / sigma[0])

plt.show()

