import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import numpy as np

n0 = 2.29 #refraction constant
lmd = 0.63e-6 #wavelength
L = 0.698 #distance to screen
l = 26e-3 #length of crystal

ringData = csvreader.readData("ring-data.csv")
k, b, dk, db = graphs.plotlsqm(ringData[0:, 0], (ringData[0:, 1] / 1000)**2, np.zeros(np.shape(ringData)[0]), 2 * ringData[0:, 1] / 1000 * 1e-3, bflag = 0) #r^2 = k * m

n0ne = lmd / l * (n0 * L)**2 / k #equasion (2)
dn0ne = lmd / l * (n0 * L)**2 / k**2 * dk #error by differential

print(f"n0 - ne = {n0ne:.4f} +- {dn0ne:.4f}")

Uhalf1 = 405
dUhalf1 = 1

print(f"Uλ/2 = ({Uhalf1} +- {dUhalf1})V by measuring max light intensity") #Uλ/2 from step 5

Uhalf2 = 960 - 480
dUhalf2 = 2

print(f"Uλ/2 = ({Uhalf2} +- {dUhalf2})V from difference in min and max signal on oscilloscope")
