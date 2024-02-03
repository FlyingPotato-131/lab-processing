import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations
import numpy as np

data = csvreader.readData('data.csv', 0)
# print(data)
Prelation = data[1, 0:] / data[2, 0:]
print(Prelation)
print()
Urelation = np.sqrt(Prelation)
print(Urelation)
print()
Ucore = 1.32 * np.sqrt(data[2, 0:])
print(Ucore)
print()
Re = 0.4021e5 * np.sqrt(data[1, 0:])
print(Re)
