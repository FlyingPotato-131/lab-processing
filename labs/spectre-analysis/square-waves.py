import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations
import math
import re
import csv
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial as poly

import glob
# files = glob.glob("data/*.csv")

def getSpectreParams(file, params = "widths", alpha = 1, plPower = 2, height = 0, startMid = False):
	data = csvreader.readTable(file, 2, titleSize = 3)
	freq = int(re.search('[0-9]+', re.search('[0-9]+Hz', file).group(0)).group(0))
	pulseT = int(re.search('[0-9]+', re.search('[0-9]+us', file).group(0)).group(0))
	coeff = alpha * freq * pulseT / 15000
	peaks = calculations.findPeaks(data, coeff, 0)
	# print(peaks)

	maximum = np.argmax(peaks[0:, 1])
	# print(maximum)
	deltaArg = 1
	while(2 * peaks[maximum + deltaArg, 1] <= peaks[maximum, 1]):
		# print(np.shape(peaks)[0] - maximum <= 1)
		deltaArg = deltaArg + 1 if np.shape(peaks)[0] - maximum >= 1 else (deltaArg - 1)
		# print(deltaArg)
	delta = peaks[maximum + deltaArg, 0] - peaks[maximum, 0]
	print(delta)

	tmp = maximum + 1
	while tmp < np.shape(peaks)[0] - 1:
		if(peaks[tmp, 0] - peaks[tmp - 1, 0] < 0.85 * delta): # or (tmp < np.shape(peaks)[0] and peaks[tmp + 1, 0] - peaks[tmp, 0] < 0.85 * delta)):
			# if(peaks[tmp + 1, 0] - peaks[tmp, 0] < 0.2 * delta):
			# 	delvalue = np.argmin(peaks[tmp : tmp + 1, 1])
			# else:
			# 	delvalue = tmp	
			# peaks = np.delete(peaks, np.argmin(peaks[tmp - 1:tmp, 1]), axis = 0)
			# if(peaks[tmp, 1] < peaks[tmp - 1, 0]):
			# 	tmp -= 1
			peaks = np.delete(peaks, tmp, axis = 0)
		else:
			tmp += 1
	# if(peaks[-1, 0] - peaks[-2, 0] > 5 * delta):
	# 	peaks = np.delete(peaks, -1, axis = 0)

	maximum = np.argmax(peaks[0:, 1])
	# print(maximum)
	tmp = maximum - 1
	while tmp > 0:
		if(peaks[tmp + 1, 0] - peaks[tmp, 0] < 0.85 * delta):
			# np.append(peaks, [peaksUnf[i]], axis = 0)
			peaks = np.delete(peaks, tmp, axis = 0)
		else:
			tmp -= 1

	specEnd = 0
	for i in range(maximum, np.shape(peaks)[0]):
		# and peaks[i, 1] <= coeff / 6
		if(i == np.shape(peaks)[0] - 1 or (peaks[i + 1, 0] - peaks[i, 0] > 1.1 * (peaks[maximum + 1, 0] - peaks[maximum, 0]))):
			specEnd = i
			break

	specBegin = 0
	for i in reversed(range(maximum)):
		# print(i)
		if(i == 0 or (peaks[i, 0] - peaks[i - 1, 0] > 1.1 * (peaks[maximum + 1, 0] - peaks[maximum, 0]))):
			specBegin = i
			break

	print(specBegin, specEnd)

	if(params == "harmonics"):
		return peaks[0 : specEnd, 0:]

	peakDist = np.average(np.array([abs(peaks[i, 0] - peaks[i - 1, 0]) for i in range(specBegin + 1, specEnd)]))
	dDist = 2 * math.sqrt(np.average(np.array([(peaks[i, 0] - peaks[i - 1, 0] - peakDist)**2 for i in range(specBegin + 1, specEnd)])) / (specEnd - specBegin - 1))

	print(peakDist, dDist)

	spec0 = poly.fit(peaks[specBegin : specEnd + 1, 0], peaks[specBegin : specEnd + 1, 1], plPower)
	t = np.linspace(peaks[specBegin, 0] - 3, peaks[specEnd, 0] + 3, num = 1000)

	if(not startMid):
		rootMax = calculations.newtonMethod(spec0, height * peaks[maximum, 1], peaks[specEnd, 0], 0.01)
		rootMin = calculations.newtonMethod(spec0, height * peaks[maximum, 1], peaks[specBegin, 0], 0.01)
	else:
		rootMax = calculations.newtonMethod(spec0, height * peaks[maximum, 1], 0.5 * peaks[specBegin, 0] + 0.5 * peaks[specEnd, 0], 0.01)
		rootMin = rootMax

	# print(rootMax, rootMin)
	spectreWidth = rootMax - max(rootMin, 0) if rootMax - rootMin > peakDist else rootMax
	R2 = np.average(np.array([abs(spec0(peaks[i, 0]) - peaks[i, 1]) for i in range(specBegin, specEnd)]))
	dWidth = abs(R2 / spec0.deriv()(spectreWidth))
	print(spectreWidth, dWidth)

	fig, ax = plt.subplots()
	plt.minorticks_on()
	plt.xlabel("ν, кГц")
	plt.ylabel("a, усл. ед.")
	plt.grid(True, "major", "both", color = "#888888")
	plt.grid(True, "minor", "both", linestyle = '--')
	ax.plot(t, spec0(t), 'r')
	ax.plot(data[:, 0], data[:, 1])
	ax.plot(peaks[:, 0], peaks[:, 1], 'r.')
	# plt.close(fig)
	plt.show()
	return peakDist, dDist, spectreWidth, dWidth, rootMin, rootMax

# for i in range(10):
# 	getSpectreParams("data/sq-1000Hz/" + str((1 + i) * 20) + "us.csv")
# for i in files:
# 	print(i)
	# getSpectreParams(i)
# getSpectreParams("data/sq-1000Hz/80us.csv")
# getSpectreParams("data/sq-100us/200Hz.csv")

# constPeaks = getSpectreParams("data/sq-1000Hz/60us.csv", "harmonics")
# print(constPeaks)

# with open("results/sq-1000Hz-60us-harmonics.csv", 'w', newline = '') as results:
# 	writer = csv.writer(results)
# 	writer.writerow(['n'] + [i + 1 for i in range(np.shape(constPeaks)[0])])
# 	writer.writerow(["nu_exp, kHz"] + constPeaks[0:, 0].tolist())
# 	writer.writerow(["nu_theor, kHz"] + [i + 1 for i in range(np.shape(constPeaks)[0])])
# 	writer.writerow(["|a|_exp"] + constPeaks[0:, 1].tolist())
# 	writer.writerow(["|a / a_1|_exp"] + [i / constPeaks[0, 1] for i in constPeaks[0:, 1]])
# 	writer.writerow(["|a / a_1|_theor"] + [abs(math.sin(math.pi * i * 60 * 10**-6 / 10**-3) / i / math.sin(math.pi * 60 * 10**-6 / 10**-3)) for i in range(1, np.shape(constPeaks)[0] + 1)])

# files1000Hz = glob.glob("data/sq-1000Hz-2/*.csv")
# widths = np.empty([0, 2])
# tau = np.empty(0)
# for dataFile in files1000Hz:
# 	print(dataFile)
# 	params = getSpectreParams(dataFile)
# 	widths = np.append(widths, [[params[2], params[3]]], axis = 0)
# 	tau = np.append(tau, int(re.search('[0-9]+', re.search('[0-9]+us', dataFile).group(0)).group(0)))

# x = np.array([1 / i / 10**-6 for i in tau])
# print(graphs.plotLsqm(x, widths[0:, 0], dy = widths[0:, 1], xlabel = "1 / τ, с^-1", ylabel = 'Δν, кГц'))
# print(widths)

# files100us = glob.glob("data/sq-100us/*Hz*.csv")
# dists = np.empty([0, 2])
# nu = np.empty(0)
# for dataFile in files100us:
# 	print(dataFile)
# 	params = getSpectreParams(dataFile)
# 	dists = np.append(dists, [[params[0], params[1]]], axis = 0)
# 	nu = np.append(nu, int(re.search('[0-9]+', re.search('[0-9]+Hz', dataFile).group(0)).group(0)))

# print(graphs.plotLsqm(nu, dists[0:, 0], dy = dists[0:, 1]))
# print(dists)

# filesCg = glob.glob("data/cg*.csv")
# cgparams = np.empty([0, 6])
# for dataFile in filesCg:
# 	print(dataFile)
# 	cgparams = np.append(cgparams, [getSpectreParams(dataFile, alpha = 0.0001)], axis = 0)

# print(cgparams)
# print(cgparams[0:, 2:4])

# filesGs = glob.glob("data/gs*.csv")
# gsparams = np.empty([0, 6])
# for dataFile in filesGs:
# 	print(dataFile)
# 	gsparams = np.append(gsparams, [getSpectreParams(dataFile, alpha = 0.01, plPower = 4, height = 0.5, startMid = True)], axis = 0)
# 	print(gsparams[-1, 5] if gsparams[-1, 4] < 0 else min(gsparams[-1, 4], gsparams[-1, 5]))

# filesAmSig = glob.glob("data/am*sig.csv")
# for dataFile in filesAmSig:
# 	print(dataFile)
# 	data = csvreader.readTable(dataFile, 2, titleSize = 3)
# 	fig, ax = plt.subplots()
# 	plt.minorticks_on()
# 	plt.grid(True, "major", "both", color = "#888888")
# 	plt.grid(True, "minor", "both", linestyle = '--')
# 	ax.plot(data[:, 0], data[:, 1])
# 	plt.show()

# filesAm = glob.glob("data/am-50kHz-2kHz/*.csv")
# adivaofm = np.empty([0, 2])

# for dataFile in filesAm:
# 	print(dataFile)
# 	data = csvreader.readTable(dataFile, 2, titleSize = 3)
# 	peaks = calculations.findPeaks(data, 1, 0)
# 	print(peaks)
# 	Amain = peaks[1, 1]
# 	Aside = 0.5 * peaks[0, 1] + 0.5 * peaks[2, 1]
# 	dataFile = dataFile.replace("data/am-50kHz-2kHz/", '')
# 	# m = float(dataFile)
# 	m = float(re.search('[+-]?([0-9]*[.])?[0-9]+', re.search('[+-]?([0-9]*[.])?[0-9]+', dataFile).group(0)).group(0))
# 	adivaofm = np.append(adivaofm, [[m, Aside / Amain]], axis = 0)

# print(adivaofm)

# print(graphs.plotLsqm(adivaofm[0:, 0], adivaofm[0:, 1], xlabel = "m", ylabel = "a_бок / a_осн"))

# filesPm = glob.glob("data/pm*.csv")

# for dataFile in filesPm:
# 	print(dataFile)
# 	data = csvreader.readTable(dataFile, 2, titleSize = 3)
# 	graphs.plot(data[0:, 0], data[0:, 1], plotFmt = '-')

data = csvreader.readTable("data/flt-400000Hz-150ns.csv", 3, titleSize = 3)

fig, ax = plt.subplots()
plt.minorticks_on()
plt.xlabel("ν, кГц")
plt.ylabel("a, усл. ед.")
plt.grid(True, "major", "both", color = "#888888")
plt.grid(True, "minor", "both", linestyle = '--')
ax.plot(data[:, 0], data[:, 1])
ax.plot(data[:, 0], data[:, 2])
# plt.close(fig)
plt.show()
fltPeaks = getSpectreParams(data[0:, 0:1])
fltPeaks = getSpectreParams("data/flt-400000Hz-150ns.csv", 3, titleSize = 3, params = "harmonics")
