import numpy as np
from numpy.polynomial import Polynomial as poly

def newtonMethod(function, value, startValue, precision):
	deriv = function.deriv()
	x = startValue
	dy = 2 * precision * abs(function(x) - value)
	while(dy > precision * abs(function(x) - value)):
		k = deriv(x)
		y = function(x) - value
		b = y - k * x
		x = -b / k
		dy = abs(y - (function(x) - value))
		# print(dy)
	return x

def closestValue(array, value):
	index = 0
	isRising = (array[1] - array[0] > 0)
	if(isRising and value < array[0]):
		return index
	while(index < np.size(array)):
		print('test')
		if((isRising and value > array[index]) or (not isRising and value < array[index])):
			if((isRising and value < array[index + 1]) or (not isRising and value > array[index + 1])):
				# break
				return index
			else:
				index += 1
	return index