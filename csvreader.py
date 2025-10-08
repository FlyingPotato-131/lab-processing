import csv
import numpy as np

def readTable(fileName, width, titleSize = 1, dataType = 'f', split = ',', enc = "utf-8"):
	file = open(fileName, newline='', encoding = enc)
	try:
		reader = csv.reader(file, delimiter = split)
		data = np.empty(0)
		i = 0
		for row in reader:
			if(i >= titleSize):
				data = np.append(data, row)
			i += 1
		data = data.reshape(np.size(data) // width, width)
		# print(data)
	finally:
		file.close()
	return data.astype(dataType)

def readMisc(fileName, width, hasTitle = False):
	return readTable(fileName, width, titleSize = hasTitle)[0]

def readMiscVert(fileName, height):
	return readTable(fileName, 1, titleSize = 0).reshape(1, height)[0]

def readData(filename, titlesize = 1, datatype = 'f', split = ',', enc = "utf-8"):
	with open(filename, newline = '', encoding = enc) as file:
		reader = csv.reader(file, delimiter = split)
		return np.array(list(reader))[titlesize:, ].astype(datatype)
