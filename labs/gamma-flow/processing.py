import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, pi

ambient = 3377

dW = 0.1

#aluminum
data = csvreader.readData("al-data.csv", 0)
data[:,1] = data[:,1] - ambient

for s in range(3):
	for i in range(1, 6):
		data[6*s+i, 0] = data[6*s+i, 0] + data[6*s+i-1, 0]
		data[6*s+i, 1] = data[6*s, 1] / data[6*s+i,1]
	data[6*s,1] = 1

dataAvg = (data[:6,:] + data[6:12,:] + data[12:,:]) / 3
dataErr = 3 * np.sqrt((data[:6,1] - dataAvg[:,1])**2 + (data[6:12,1] - dataAvg[:,1])**2 + (data[12:,1] - dataAvg[:,1])**2)

k, b, dk, db = graphs.plotlsqm(dataAvg[:, 0], np.log(dataAvg[:, 1]), np.ones(6) * dW, dataErr / dataAvg[:, 1], title = "Логарифмическое ослабление потока, Al", xlabel = "толщина образца, мм", ylabel = "ln(N0 / N)", bflag = False)
print(f"μ_Al = ({k * 10:.5f} +- {dk * 10:.5f})cm^-1")
print(f"hω_Al = 0.75 MeV")
print()

#iron
data = csvreader.readData("fe-data.csv", 0)
data[:,1] = data[:,1] - ambient

for s in range(3):
	for i in range(1, 6):
		data[6*s+i, 0] = data[6*s+i, 0] + data[6*s+i-1, 0]
		data[6*s+i, 1] = data[6*s, 1] / data[6*s+i,1]
	data[6*s,1] = 1

dataAvg = (data[:6,:] + data[6:12,:] + data[12:,:]) / 3
dataErr = 3 * np.sqrt((data[:6,1] - dataAvg[:,1])**2 + (data[6:12,1] - dataAvg[:,1])**2 + (data[12:,1] - dataAvg[:,1])**2)

k, b, dk, db = graphs.plotlsqm(dataAvg[:, 0], np.log(dataAvg[:, 1]), np.ones(6) * dW, dataErr / dataAvg[:, 1], title = "Логарифмическое ослабление потока, Fe", xlabel = "толщина образца, мм", ylabel = "ln(N0 / N)", bflag = False)
print(f"μ_Fe = ({k * 10:.4f} +- {dk * 10:.4f})cm^-1")
print(f"hω_Fe = 0.77 MeV")
print()

#lead
data = csvreader.readData("pb-data.csv", 0)
data[:,1] = data[:,1] - ambient

for s in range(3):
	for i in range(1, 5):
		data[5*s+i, 0] = data[5*s+i, 0] + data[5*s+i-1, 0]
		data[5*s+i, 1] = data[5*s, 1] / data[5*s+i,1]
	data[5*s,1] = 1

dataAvg = (data[:5,:] + data[5:10,:] + data[10:,:]) / 3
dataErr = 3 * np.sqrt((data[:5,1] - dataAvg[:,1])**2 + (data[5:10,1] - dataAvg[:,1])**2 + (data[10:,1] - dataAvg[:,1])**2)

k, b, dk, db = graphs.plotlsqm(dataAvg[:, 0], np.log(dataAvg[:, 1]), np.ones(5) * dW, dataErr / dataAvg[:, 1], title = "Логарифмическое ослабление потока, Pb", xlabel = "толщина образца, мм", ylabel = "ln(N0 / N)", bflag = False)
print(f"μ_Pb = ({k * 10:.3f} +- {dk * 10:.3f})cm^-1")
print(f"hω_Pb = 0.75 MeV")
print()

avgE = (0.75 + 0.75 + 0.77) / 3
errE = sqrt((0.75 - avgE)**2 + (0.77 - avgE)**2 + (0.75 - avgE)**2)

print(f"hω_avg = ({avgE:.3f} +- {errE:.3f}) MeV")
