import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations
import math
import numpy as np

def sumall(*args):
    s = 0

    for num in args:
        s += float(num)

    return s

def rnderror(*args):
    sqrsum = 0
    ln = len(args)
    avg = sumall(*args)/ln

    for num in args:
        sqrsum += (float(num) - avg)**2

    return math.sqrt(sqrsum/ln/(ln-1))

def avgerror(*args):
    return sumall(*args)/len(args), rnderror(*args)

data = csvreader.readData("data.csv")
T, dT0 = avgerror(*data[:, 1])
dT = math.sqrt(dT0**2 + 6**2)
print(f"T = ({T:.0f} +- {dT:.0f}) K")
print(1 / (1 / data[:, 1] + 0.65e-6 / 1.44e-2 * np.log(0.374 * math.sqrt(51e-9) * data[:, 1])))
