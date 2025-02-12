import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as sci
from math import sqrt, pi

#constant definitions
N = 50 #number of points
gamma = 7/5 #assuming air as two atom molecule
F2 = 1 #cross section of high pressure area
F1 = 5 #cross section of low pressure area
M5 = 1 #mach number at nozzle entrance
Fs = np.linspace(1, 6.5, N) #cross section at shockwave in nozzle (parameter)
k = 3655 #pressure sensor coefficient
rho = 13596 #mercury density
g0 = 9.81 #gravity acceleration

#function definitions
def q(M): #cross section fraction depending on mach (q function, page 27)
	return M * (2 / (gamma + 1) * (1 + (gamma - 1) / 2 * M**2))**(-(gamma + 1) / 2 / (gamma - 1))

def tau(M): #~temperature on shockwave (tau function)
	return 1 / (1 + (gamma - 1) / 2 * M**2)

def pi(M): #~pressure on shockwave (pi function)
	return (1 + (gamma - 1) / 2 * M**2)**(-gamma / (gamma - 1))

def M_eq(M, M_, F_, F): #define equation for finding M on lavalle nozzle, q(M) / q(M_) = F_ / F
	return q(M) - F_ / F * q(M_)

def g(pp): #=u2/a1, equation 8.3
	return (pp - 1) / np.sqrt((gamma * (gamma - 1) / 2) * (1 + (gamma + 1) / (gamma - 1) * pp))

def p21_eq(pp, h): #define equation for finding p2/p1, g(p2 / p1) = h(Fs)
	return g(pp) - h

#calculate mach numbers
M6 = sci.fsolve(M_eq, 2 * np.ones(N), args=(M5, F2, Fs)) #solve q function equation

M7 = (2 + (gamma - 1) * M6) / (2 * gamma * M6 - (gamma - 1)) #shockwave mach number relation

M3 = sci.fsolve(M_eq, 0.5 * np.ones(N), args=(M7, Fs, F1)) #solve q function equation

#calculate speed of sound fractions
#a4 / a1 = 1 since T4 = T1
a37 = np.sqrt(tau(M3) / tau(M7)) #a3 / a7 from tau functions <=> temperature on shockwave
a76 = np.sqrt(tau(M7) / tau(M6)) #a7 / a6, same method
a65 = np.sqrt(tau(M6) / tau(M5)) #a6 / a5, same method
a54 = 2 / (gamma + 1) #page 28

hFs = M3 * a37 * a76 * a65 * a54 #h(Fs) = u2 / a1, g(p2/p1) = h(Fs)

#calculate pressure fractions
p37 = pi(M3) / pi(M7) #p3 / p7 from pi functions
p76 = (1 + gamma * M6**2) / (1 + gamma * M7**2) #p7 / p6 from pressure on shockwave
p65 = pi(M6) / pi(M5) #p6 / p5 from pi functions
p54 = 0.279 #p5 / p4 from adiabatic flow

fFs = p37 * p76 * p65 * p54 #f(Fs), p2/p1 = p4/p1 * f(Fs)

#solve for p2/p1(p4/p1)
p21 = sci.fsolve(p21_eq, np.ones(N), args=(hFs)) #p2/p1
p41 = p21 / fFs #p4/p1

data = csvreader.readData("data.csv")
p41_exp = (data[:, 1] / 250 * 6 + 1) / (data[:, 0]*1e-3 * rho * g0 / 101330) #p4/p1 measured
dp41_exp = (1 / 250 * 6) / (data[:, 0]*1e-3 * rho * g0 / 101330) + (data[:, 1] / 250 * 6 + 1) / (data[:, 0]*1e-3 * rho * g0 / 101330)**2 * (1e-3 * rho * g0 / 101330) #p4/p1 error

pdiff = data[:, 2] / k * 1e6 #p2 - p1 from pressure sensor, Pa
dpdiff = 0.1 / k * 1e6 #p2 - p1 error
p21_se = (data[:, 0]*1e-3 * rho * g0 + pdiff) / (data[:, 0]*1e-3 * rho * g0) #p2/p1 from pressure sensor
dp21_se = (1e-3 * rho * g0 + dpdiff) / (data[:, 0]*1e-3 * rho * g0) + (data[:, 0]*1e-3 * rho * g0 + pdiff) / (data[:, 0]*1e-3 * rho * g0)**2 * (1e-3 * rho * g0) #p2/p1 error

Vs = 0.05 / (data[:, 4] - data[:, 3]) * 1e6 #shockwave velocity
dVs = 0.05 / (data[:, 4] - data[:, 3])**2 * 1e6 * 2 #shockwave velocity error
Ms = Vs / np.sqrt(gamma * 8.31 * (24 + 273) / 29e-3) #shockwave mach
dMs = dVs / np.sqrt(gamma * 8.31 * (24 + 273) / 29e-3) #shockwave mach error
p21_ls = (Ms**2 - (gamma - 1) / 2 / gamma) / (gamma + 1) * 2 * gamma #p2/p1 from laser sensor
dp21_ls = (2 * Ms * dMs) / (gamma + 1) * 2 * gamma #p2/p1 from laser sensor

fig, ax = graphs.basePlot()

ax.plot(p41, p21, label = "theor", linewidth=1)
ax.errorbar(p41_exp, p21_se, dp21_se, dp41_exp, fmt = '.', label = "pressure sensor", linewidth=1)
ax.errorbar(p41_exp, p21_ls, dp21_ls, dp41_exp, fmt = '.', label = "laser sensor", linewidth=1)

plt.title("p2/p1 (p4/p1)")
plt.xlabel("p4/p1")
plt.ylabel("p2/p1")

plt.legend()

plt.show()
