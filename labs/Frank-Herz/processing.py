import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

from math import pi
import math
import numpy as np
import matplotlib.pyplot as plt

data4 = csvreader.readData("data4.csv", 0)
data5 = csvreader.readData("data5.csv", 0)
data8 = csvreader.readData("data8.csv", 0)

fig, ax = graphs.basePlot()
plt.title("ВАХ")
plt.xlabel("V_a, V")
plt.ylabel("I_k, mA")
ax.errorbar(data4[:, 0], data4[:, 1], 0.005, 0.01, label = "4v")
ax.errorbar(data5[:, 0], data5[:, 1], 0.005, 0.01, label = "5v")
ax.errorbar(data8[:, 0], data8[:, 1], 0.005, 0.01, label = "8v")
plt.legend()

print(f"ΔV = ({15.0} +- {1.6}) V")

plt.show()