import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import numpy as np
import matplotlib.pyplot as plt

h = 30e-3
w = 6e-3
I = w * h**3 / 12

def floatwith0(string):
	if(string == ''):
		return 0
	else:
		return float(string)

raw = csvreader.readData('data.csv', datatype = 'str')
data = np.empty(np.shape(raw))

for h in range(np.shape(raw)[0]):
	for w in range(np.shape(raw)[1]):
		data[h, w] = floatwith0(raw[h, w])

# print(data)
for i in range(9, 13):
	partial = data[i, 6:20]
	# print(partial)
	dist = -0.5 * (partial[::3] + partial[1::3])
	tension = h * data[i, 0] * 11 * 0.1 * 9.8 / 4 / I * dist
	fig, ax = graphs.basePlot()
	k, b, dk, db = graphs.lsqm(np.arange(2, 7), tension, np.zeros(np.size(tension)), partial[2::3], bflag = 1)

	print(f'P = {data[i, 0] * 11 * 0.1 * 9.8:.2f} n')
	print(f'k = {k * 1e-9:.1f} +- {dk * 1e-9:.1f} * 10^9 Pa')
	print(f'b = {b * 1e-9:.1f} +- {db * 1e-9:.1f} * 10^9 Pa\n')

	plt.title(f'σ of i, P = {data[i, 0] * 11 * 0.1 * 9.8:.2f} n')
	plt.xlabel('i')
	plt.ylabel('σ, Pa')

	ax.plot(np.arange(2, 7), tension, '.')
	ax.plot(np.arange(2, 7), k * np.arange(2, 7) + b)
	plt.show()

# for i in range(1, 4):
# 	# partial = data[data[0:, i] != 0]
# 	# partial = data[np.all(data != 0, axis = 1)]
# 	partial = np.array([row for row in data if row[3 + 3 * i] != 0])
# 	# print(partial, '\n')
# 	dist = (0.5 * partial[0:, 3 + 3 * i] + 0.5 * partial[0:, 3 + 3 * i + 1]) * 1e-3
# 	tension = h * partial[0:, 0] * 11 * 0.1 / 4 / I * dist
# 	fig, ax = graphs.basePlot()

# 	k, b, dk, db = graphs.lsqm(partial[0:, 0], tension, np.zeros(np.size(tension)), partial[0:, 3 + 3 * i + 2], bflag = 1)
# 	print(f'k = {k} +- {dk}')
# 	print(f'b = {b} +- {db}\n')

# 	plt.title('σ of i')
# 	plt.xlabel('i')
# 	plt.ylabel('σ, Pa')

# 	ax.plot(partial[0:, 0], tension, '.')
# 	ax.plot(partial[0:, 0], k * partial[0:, 0] + b)
# 	plt.show()
