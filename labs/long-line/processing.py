import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import math
from math import pi
import numpy as np

resonancedata = csvreader.readData('resonance-data.csv')

L = 5030
d = 0.137
harmonics = np.array([1, 2, 3, 4])
vpharray = np.empty(0)
dvpharray = np.empty(0)

for i in range(4):
	k, b, dk, db = graphs.plotlsqm(harmonics, L * resonancedata[0:, i], np.zeros(4), L * 0.01e6 * np.ones(4) + 0.1 * resonancedata[0:, i] * np.ones(4), bflag = False, title = 'Lν of harmonic')
	vpharray = np.append(vpharray, k)
	dvpharray = np.append(dvpharray, dk)

vph, dvph = calculations.avgerror(vpharray, dvpharray)
print(f'Vф = {vph * 1e-6:.0f} +- {dvph * 1e-6:.0f} * 10^6 cм/с')

phasedata = csvreader.readData('phase-data.csv')
alphaw = 1 / L * np.log(phasedata[0:, 3] / phasedata[0:, 4])
dalphaw = np.sqrt((1 / L**2 * np.log(phasedata[0:, 3] / phasedata[0:, 4]) * 0.1)**2 + (1 / L / phasedata[0:, 3] * 0.01)**2 + (1 / L / phasedata[0:, 4] * 0.01)**2)
phaseshift = 2 * pi * phasedata[0:, 1] * phasedata[0:, 0] + phasedata[0:, 2] * pi
dphaseshift = np.sqrt((2 * pi * 0.01e-7 * phasedata[0:, 0])**2 + (2 * pi * phasedata[0:, 1] * 0.01e6)**2)
kw = phaseshift / L
dkw = dphaseshift / L

x1 = (2 * pi * phasedata[0:, 0])**2
dx1 = 4 * pi**2 * phasedata[0:, 0] * 0.01e6
y1 = kw**2 - alphaw**2
dy1 = kw * dkw + alphaw * dalphaw
graphs.plotlsqm(x1 * 1e-15, y1, dx1 * 1e-15, dy1, bflag = False, title = 'k^2 - α^2 of ω^2', xlabel = 'ω^2, Hz^2 * 10^15', ylabel = 'k^2 - α^2, cm^-2')
k, b, dk, db = graphs.lsqm(x1, y1, dx1, dy1, bflag = False)

LxCx = k * 3e10**2 
dLxCx = dk * 3e10**2

R0 = 7.23e-11 #в сгс
LxdivCx = (R0 * 3e10)**2

Lx = math.sqrt(LxCx * LxdivCx)
dLx = 0.5 * math.sqrt(LxdivCx / LxCx) * dLxCx
print(f'Lx = {Lx:.2f} +- {dLx:.2f} сгс')
Cx = math.sqrt(LxCx / LxdivCx)
dCx = 0.5 * math.sqrt(1 / LxdivCx / LxCx) * dLxCx
print(f'Cx = {Cx:.3f} +- {dCx:.3f} сгс')

x2 = np.sqrt(phasedata[0:, 0])
dx2 = 1 / np.sqrt(phasedata[0:, 0]) * 0.01e6
y2 = alphaw
dy2 = dalphaw
k, b, dk, db = graphs.plotlsqm(x2, y2, dx2, dy2, bflag = False, title = 'α of sqrt(ν)', xlabel = 'sqrt(ν), sqrt(Hz)', ylabel = 'α, cm^-1')
sigma = (2 / math.sqrt(LxdivCx) / d / k)**2
dsigma = (2 / math.sqrt(LxdivCx) / d / k) * 2 / math.sqrt(LxdivCx) / d / k**2 * dk
# print(sigma, dsigma)
print(f'σ_A = {sigma / 1e15:.0f} +- {dsigma / 1e15:.0f} * 10^15 сгс')

x3 = (np.sqrt(phasedata[0:, 0]))**3
y3 = alphaw * kw
k, b, dk, db = graphs.plotlsqm(x3, y3, bflag = False, title = 'kα of ν^3/2', xlabel = 'ν^3/2, Hz^3/2', ylabel = 'kα, cm^-2')
sigma = (4 * pi * Cx / d / 3e10 / k)**2
dsigma = (4 * pi * Cx / d / 3e10 / k) * (4 * pi * dCx / d / 3e10 / k + 4 * pi * Cx / d / 3e10 / k**2 * dk)
print(f'σ_Б = {sigma / 1e15:.0f} +- {dsigma / 1e15:.0f} * 10^15 сгс')
