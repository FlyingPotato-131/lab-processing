import csv
import numpy as np

def readTable(fileName, width, titleSize = 1, dataType = 'f'):
	file = open(fileName, newline='')
	try:
		reader = csv.reader(file)
		data = np.empty(0)
		for row in reader:
			data = np.append(data, row)
		data = data.reshape(np.size(data) // width, width)
		# print(data)
	finally:
		file.close()
	return data[titleSize:, ].astype(dataType)

def readMisc(fileName, width, hasTitle = False):
	return readTable(fileName, width, titleSize = hasTitle)[0]

def readMiscVert(fileName, height):
	return readTable(fileName, 1, titleSize = 0).reshape(1, height)[0]

def readData(filename, titlesize = 1, datatype = 'f'):
	with open(filename, newline = '') as file:
		reader = csv.reader(file)
		return np.array(list(reader))[titlesize:, ].astype(datatype)
