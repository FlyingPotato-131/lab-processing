import math
import numpy as np
from numpy.polynomial import Polynomial as poly

def newtonMethod(function, value, startValue, precision, rootMin = None, rootMax = None):
	deriv = function.deriv()
	x = startValue
	dy = 2 * precision * abs(function(x) - value)
	while(dy > precision * abs(function(startValue) - value)):
		k = deriv(x)
		y = function(x) - value
		b = y - k * x
		x = -b / k
		dy = abs(y - (function(x) - value))
		# if(dy <= precision * abs(function(startValue) - value) and not(rootMin is None) and x < rootMin):
		# 	return rootMin
		# if(dy <= precision * abs(function(startValue) - value) and not (rootMax is None) and x > rootMax):
		# 	return rootMax
	return min(rootMax, max(rootMin, x)) if rootMin != None and rootMax != None else x
	# return x

def closestValue(array, value):
	index = 0
	isRising = (array[1] - array[0] > 0)
	if((isRising and value < array[0]) or (not isRising and value > array[0])):
		return index
	while(index < np.size(array)):
		# print('test')
		# print(isRising)
		# print(index)
		if((isRising and value > array[index]) or (not isRising and value < array[index])):
			if((isRising and value < array[index + 1]) or (not isRising and value > array[index + 1])):
				# break
				return index
			else:
				index += 1
	return index

def findPeaks(data, coeff, lowerBound):
	peaks = np.empty([0, 2])
	for i in range(1, np.shape(data)[0] - 1):
		if((data[i, 1] - data[i-1, 1]) / (data[i, 0] - data[i-1, 0]) > coeff and data[i+1, 1] < data[i, 1] and data[i, 1] >= lowerBound):
			peaks = np.append(peaks, [data[i]], axis = 0)
		elif((data[i, 1] - data[i+1, 1]) / (data[i+1, 0] - data[i, 0]) > coeff and data[i-1, 1] < data[i, 1] and data[i, 1] >= lowerBound):
			peaks = np.append(peaks, [data[i]], axis = 0)
	return peaks

def R2(function, x, y):
	return math.sqrt(np.sum(np.array([(y[i] - function(x[i]))**2 for i in range(len(x))])))

def avgerror(array, darray):
	avg = np.average(array)
	return [avg, math.sqrt(np.average((array - avg)**2) + np.average(darray)**2)]