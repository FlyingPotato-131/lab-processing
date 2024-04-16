import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations
import numpy as np
import matplotlib.pyplot as plt
import math

besselData = csvreader.readData("bessel-data.csv")
for i in range(np.shape(besselData)[0]):
	print(f"f_{i+1} = ({(besselData[i, 5]**2 - besselData[i, 3]**2) / 4 / besselData[i, 5]:.2f} +- {math.sqrt(((1 / 4 - besselData[i, 3]**2 / 4 / besselData[i, 5]**2) * 1)**2 + (besselData[i, 3] / 2 / besselData[i, 5] * 1)**2):.2f}) mm")
print()

abbeData = csvreader.readData("abbe-data.csv")
for i in range(np.shape(abbeData)[0]):
	print(f"f_{i+1} = ({abbeData[i, 3] / (abbeData[i, 0] / abbeData[i, 4] - abbeData[i, 1] / abbeData[i, 4]):.0f} +- {np.sqrt((1 / (abbeData[i, 0] / abbeData[i, 4] - abbeData[i, 1] / abbeData[i, 4]))**2 + (abbeData[i, 3] / (abbeData[i, 0] - abbeData[i, 1]))**2 + (abbeData[i, 3] / (abbeData[i, 0] / abbeData[i, 4] - abbeData[i, 1] / abbeData[i, 4])**2 / abbeData[i, 4])**2):.0f}) mm")
