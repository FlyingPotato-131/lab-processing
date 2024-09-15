import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations
import numpy as np

data1 = np.array([2, 5, 6, 8, 9, 10, 11, 12])
print(data1 * 45 / 1085)
x = np.array([1, 2, 3, 4, 5, 6, 7, 8]) - 1
print(graphs.plotlsqm(x, (data1 * 45 / 1085)**2, np.zeros(8), 2 * data1 * 45 / 1085 * np.ones(8) * 45 / 1085, bflag = 0, title = "r^2 over m", xlabel = "m", ylabel = "r^2, mm^2"))