import numpy as np
import matplotlib.pyplot as plt
import math
from numpy.polynomial import Polynomial as poly

def baseLsqm(x, y):
	if(np.size(x) != np.size(y)):
		raise ValueError('size of axes not identical')
	x2 = np.array([n**2 for n in x])
	y2 = np.array([n**2 for n in y])
	xy = np.array([x[i] * y[i] for i in range(np.size(x))])
	xavg = np.average(x)
	yavg = np.average(y)
	x2avg = np.average(x2)
	y2avg = np.average(y2)
	xyavg = np.average(xy)
	k = (xyavg - xavg * yavg) / (x2avg - xavg ** 2)
	b = yavg - k * xavg
	dk = math.sqrt(((y2avg - yavg ** 2) / (x2avg - xavg**2) - k**2) / np.size(x))
	db = dk * math.sqrt(x2avg - xavg**2)
	return [k, b, dk, db]

def lsqm(x, y, dx = None, dy = None):
	if(dx.any() == None and dy.any() == None):
		return baseLsqm(x, y)
	if(dx.any() == None):
		dx = np.zeros(np.size(x))
	if(dy.any() == None):
		dy = np.zeros(np.size(y))
	k, b, dk, db = baseLsqm(x, y)
	Dk = math.sqrt(dk**2 + (dy[-1] / x[-1])**2 + ((y[-1] - b) * dx[-1] / x[-1]**2)**2)	
	Db = math.sqrt(db**2 + dy[-1]**2 + (k * dx[-1])**2)
	return [k, b, Dk, Db]

def plot(x, y, dx = None, dy = None, filename = None, plotFmt = '.', title = None, xlabel = None, ylabel = None):
	fig, ax = plt.subplots()
	plt.minorticks_on()
	plt.grid(True, "major", "both", color = "#888888")
	plt.grid(True, "minor", "both", linestyle = '--')
	plt.title(title)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)	
	ax.errorbar(x, y, dx, dy, fmt = plotFmt)
	if(filename != None):
		plt.savefig(filename)
	plt.show()

def plotLsqm(x, y, dx = None, dy = None, filename = None, plotFmt = '.', lsqmFmt = 'r', title = None, xlabel = None, ylabel = None):
	k, b, dk, db = lsqm(x, y, dx, dy)
	fig, ax = plt.subplots()
	plt.minorticks_on()
	plt.grid(True, "major", "both", color = "#888888")
	plt.grid(True, "minor", "both", linestyle = '--')
	plt.title(title)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)	
	ax.plot([x[0], x[-1]], [k * x[0] + b, k * x[-1] + b], lsqmFmt)
	ax.errorbar(x, y, dx, dy, fmt = plotFmt)	
	if(filename != None):
		plt.savefig(filename)
	plt.show()
	return [k, b, dk, db]

def plotPoly(n, x, y, dx = None, dy = None, filename = None, plotFmt = '.', polyFmt = 'r', title = None, xlabel = None, ylabel = None):
	fig, ax = plt.subplots()
	plt.minorticks_on()
	plt.grid(True, "major", "both", color = "#888888")
	plt.grid(True, "minor", "both", linestyle = '--')
	plt.title(title)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)	
	approx = poly.fit(x, y, n)
	t = np.linspace(x[0], x[-1], num = 1000)
	ax.plot(t, approx(t), polyFmt)
	ax.errorbar(x, y, dy, dx, fmt = plotFmt)
	if(filename != None):
		plt.savefig(filename)
	plt.show()
	return approx