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

params = csvreader.readTable('params.csv', 11)

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
	Hmax = params[i, 7]/2 * params[i, 5] * params[i, 0] / 0.02 / params[i, 3]
	DHmax = 1 * params[i, 5] * params[i, 0] / 0.2 / params[i, 3]

	Hc = params[i, 9]/2 * params[i, 5] * params[i, 0] / 0.02 / params[i, 3]
	DHc = 1 * params[i, 5] * params[i, 0] / 0.2 / params[i, 3]

	Bmax = params[i, 8]/2 * tau / sn * params[i, 6]
	DBmax = 1 * tau / sn * params[i, 6]

	Br = params[i, 10]/2 * tau / sn * params[i, 6]
	DBr = 1 * tau / sn * params[i, 6]

	muMax = muMaxRaw * params[i, 6] * tau * R0 * params[i, 3] / sn / params[i, 5] / params[i, 0]
	DmuMax = DmuMaxRaw * params[i, 6] * tau * R0 * params[i, 3] / sn / params[i, 5] / params[i, 0]

	mu0 = mu0Raw * params[i, 6] * tau * R0 * params[i, 3] / sn / params[i, 5] / params[i, 0]
	Dmu0 = Dmu0Raw * params[i, 6] * tau * R0 * params[i, 3] / sn / params[i, 5] / params[i, 0]

	print('mu_max', muMax, DmuMax)
	print('mu_0', mu0, Dmu0)
	print()
	print('Hmax', Hmax, DHmax)	
	print('Hc', Hc, DHc)
	print('Bmax', Bmax, DBmax)
	print('Br', Br, DBr)

	ax.plot(lim1(y), y)
	ax.plot(lim1(Bm), Bm, 'o')
	ax.plot(data[0:, 1], data[0:, 0], '.')
	
	plt.show()

