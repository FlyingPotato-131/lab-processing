import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations
import math
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import Polynomial as poly
# import scipy

tau = 400
R0 = 0.2

params = csvreader.readTable('params.csv', 7)

for i in range(3):
	sn = params[i - 1, 2] * params[i - 1, 1]

	data = csvreader.readTable('data' + str(i+1) + '.csv', 2)

	fig, ax = plt.subplots()

	plt.title('Начальная кривая намагничивания')
	plt.xlabel('H, деления на осциллографе')
	plt.ylabel('B, деления на осциллографе')

	plt.minorticks_on()
	plt.grid(True, "major", "both", color = "#888888")
	plt.grid(True, "minor", "both", linestyle = '--')

	y = np.linspace(data[0, 0], 0 * data[-1, 0], num = 1000)
	lim1 = poly.fit(data[0:, 0], data[0:, 1], 3 if i == 2 else 5)

	# R2 = calculations.R2(lim1, data[0:, 0], data[0:, 1])

	D1 = lim1.deriv()

	R2 = np.average(np.array([(D1(0.5 * data[i, 0] + 0.5 * data[i+1, 0]) - (data[i+1, 1] - data[i, 1]) / (data[i+1, 0] - data[i, 0])) if data[i+1, 0] != data[i, 0] else 0 for i in range(np.size(data[0:, 0]) - 1)]))
	print(R2)
	print()

	D2 = D1.deriv()
	Bm = calculations.newtonMethod(D2, 0, -0 + 0 * data[0, 1] + 0 * data[-1, 1], 0.0001, data[-1, 1], data[0, 1])

	muMaxRaw = 1 / D1(Bm)
	DmuMaxRaw = 1 / D1(Bm) / D1(0.5 * data[0, 0] + 0.5 * data[-1, 0]) * R2
	# R2 / D2(Bm)
	mu0Raw = 1 / D1(0)
	Dmu0Raw = 1 / D1(0) / D1(0.5 * data[0, 0] + 0.5 * data[-1, 0]) * R2
	# print(muMaxRaw)
	muMax = muMaxRaw * params[i-1, 6] * tau * R0 * params[i-1, 3] / sn / params[i-1, 5] / params[i-1, 0]
	DmuMax = DmuMaxRaw * params[i-1, 6] * tau * R0 * params[i-1, 3] / sn / params[i-1, 5] / params[i-1, 0]
	mu0 = mu0Raw * params[i-1, 6] * tau * R0 * params[i-1, 3] / sn / params[i-1, 5] / params[i-1, 0]
	Dmu0 = Dmu0Raw * params[i-1, 6] * tau * R0 * params[i-1, 3] / sn / params[i-1, 5] / params[i-1, 0]

	print(muMax, DmuMax)
	print(mu0, Dmu0)
	print()

	ax.plot(lim1(y), y)
	ax.plot(lim1(Bm), Bm, 'o')
	ax.plot(data[0:, 1], data[0:, 0], '.')
	
	plt.show()
