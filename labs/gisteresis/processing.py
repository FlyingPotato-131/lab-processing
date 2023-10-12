import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations
import math
import matplotlib.pyplot as plt
import numpy as np
# from numpy.polynomial import Polynomial as poly
import scipy

data = csvreader.readTable("data.csv", 6, titleSize = 3)
def arrayexp(array):
	# return map(math.exp, array)
	return np.array([math.exp(x) for x in array])

def arrayatan(array):
	return np.array([math.atan(x) for x in array])

def fapprox(x, m, a, b, c, n):
	# return a * math.exp(-b * x ** (-c))
	# return a * math.e**(-b * x**(-c))
	# return a * arrayexp(-(abs(x / b))**(-c)) + d * x
	return m * (arrayatan(a * x + b) + c)**n

# print(data)

for i in range(3):
	fig, ax = plt.subplots()

	plt.minorticks_on()
	plt.grid(True, "major", "both", color = "#888888")
	plt.grid(True, "minor", "both", linestyle = '--')

	x = np.linspace(data[0, 2 * i + 1], 0 * data[-1, 2 * i + 1], num = 1000)
	lim1 = scipy.optimize.curve_fit(fapprox, data[0:, 2 * i + 1], data[0:, 2 * i])
	# print(lim1)
	ax.plot(x, fapprox(x, *(lim1[0])))
	# lim1 = poly.fit(data[0:, 2 * i], data[0:, 2 * i + 1], 3)
	# ax.plot(lim1(y1), y1)
	# x1 = np.linspace(data[0, 2 * i + 1] + 1, 0 * data[-1, 2 * i + 1], num = 1000)
	# lim1 = poly.fit(data[0:, 2 * i + 1], data[0:, 2 * i], 2)
	# D1 = np.deriv(lim1)
	# D2 = np.deriv(D1)
	# D1 = lim1.deriv()
	# D2 = D1.deriv()
	# mu0 = D1(data[0, 2 * i + 1])
	# muM = calculations.newtonMethod(D2, 0, -0 + 0 * data[0, 2 * i + 1] + 0 * data[-1, 2 * i + 1], 0.0001, data[-1, 2 * i + 1], data[0, 2 * i + 1])

	# ax.plot(x, fapprox(x, lim1))
	# ax.plot(lim1(y1), y1)
	# ax.plot(x1, D1(x1))
	# ax.plot(x1, D2(x1))

	ax.plot(data[0:, 2 * i + 1], data[0:, 2 * i], '.')

	# ax.plot(data, lim1(mu0), 'o')
	# ax.plot(H_max, lim1(H_max), 'o')
	# ax.plot(lim1(muM), muM, 'o')
	
	plt.show()
