import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import math
import matplotlib.pyplot as plt
import numpy as np
import scipy

p = 2 * math.sin(math.radians(10) / 2) #constant in eq 7.4
l0 = 632.8e-6 #wavelength in mm
v0 = 6.06 #mm/s
# v0 = 4 #mm/s

data = csvreader.readData("data.csv", titlesize = 0)
samplesize = [0, 8, 8, 9, 9, 8, 7] #amount of points for each x

delta = np.empty(6)

for i in range(1, 7): #plot separate velocity profiles
	fig, ax = graphs.basePlot()
	vaxis = np.linspace(data[1, 0] * l0 / p, data[samplesize[i], 0] * l0 / p)

	ax.errorbar(data[1:samplesize[i]+1, 0] * l0 / p, data[1:samplesize[i]+1, i] * 0.01, 0.01, 100 * l0 / p, '.', label = "measured profile")
	curve = np.poly1d(np.polyfit(data[1:samplesize[i]+1, 0] * l0 / p, data[1:samplesize[i]+1, i] * 0.01, 4))
	delta[i-1] = curve(v0)
	ax.plot(vaxis, curve(vaxis))


	ax.plot([v0, v0], [data[1, i] * 0.01, data[samplesize[i], i] * 0.01], '--', label = "average flow velocity")
	plt.xlabel("v, mm/s")
	plt.ylabel("y, mm")
	plt.title(f"velocity profile, x = {data[0, i] - 39.573212315627956}")
	plt.legend()
	# k, b, dk, db = graphs.lsqm()	
	plt.show()

fig, ax = graphs.basePlot()
xaxis = np.linspace(60, 84, 100)
ax.plot(data[0, 1:], delta, '.', label = "experimental")
k, b, dk, db = graphs.lsqm(data[0, 1:], delta**2)
ax.plot(xaxis, np.sqrt(k * xaxis + b), label = "approximation")
plt.xlabel("x, mm")
plt.ylabel("delta, mm")
plt.title("boundary layer profile")
plt.legend()
plt.show()

x0 = -b / k
print(x0)

def eq_rhs(eta, y):
	d0, d1, d2 = y
	return np.array([d1, d2, -0.5 * d0 * d2])

fig, ax = graphs.basePlot()
for i in range(1, 6): #plot velocity profiles with theoretical curve
	param = 0.5
	dparam = 0.25 #parameter change on iteration
	error = scipy.integrate.solve_ivp(eq_rhs, [0, 10], [0, 0, param]).y[1, -1] - 1
	while abs(error) > 0.001:
		if error > 0:
			param = param - dparam
		else:
			param = param + dparam
		dparam = dparam / 2
		error = scipy.integrate.solve_ivp(eq_rhs, [0, 10], [0, 0, param]).y[1, -1] - 1
	solution = scipy.integrate.solve_ivp(eq_rhs, [0, 10], [0, 0, param], max_step = 0.01)
	ax.scatter(data[1:samplesize[i]+1, i] * 0.01 / math.sqrt(data[0, i] - x0), data[1:samplesize[i]+1, 0] * l0 / p / v0, label = f"x = {data[0, i]}") #TODO: figure out why this is so fucked
	ax.plot(solution.t, solution.y[1, :])
plt.xlabel("eta")
plt.ylabel("v/v0")
plt.title("undimensioned velocity profiles")
plt.show()
