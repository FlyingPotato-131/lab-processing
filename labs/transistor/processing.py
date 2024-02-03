import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations
import numpy as np
import matplotlib.pyplot as plt

R1 = 0.3
R2 = 0.68
beta = np.empty(0)
freq = np.empty(0)
dbeta = np.empty(0)

data = csvreader.readData('data.csv', titlesize = 0)
datavb = csvreader.readData('data-vb.csv', titlesize = 0)

for i in range(1, np.size(data[0, 0:])):
	f = data[0, i]
	k, b, dk, db = graphs.plotlsqm(abs(data[1:, 0] - datavb[1:, i]), data[1:, i], np.ones(np.size(data[1:, 0])), np.ones(np.size(data[1:, i])), title = f'Uout of Uin at {f} kHz', xlabel = 'Uin, mV', ylabel = 'Uout, mV')
	beta = np.append(beta, k * R1 / R2)
	dbeta = np.append(dbeta, dk * R1 / R2)
	freq = np.append(freq, f)

print(*[f'{B:.2f}' for B in beta], sep = ',')
print(*[f'{B:.2f}' for B in dbeta], sep = ',')

fig, ax = graphs.basePlot()
plt.title('beta of freq')
plt.xlabel('freq, kHz')
plt.ylabel('beta')
k, b, dk, db = graphs.lsqm(freq[3:], beta[3:])
ax.plot(freq[3:], k * freq[3:] + b)
ax.errorbar(freq, beta, dbeta, 0, '.')
plt.show()

print(f'k1 = {k:.4f} +- {dk:.4f} ms')

# fig, ax = graphs.basePlot()
# plt.title('ln beta of ln freq')
# plt.xlabel('ln freq')
# plt.ylabel('ln beta')
# k, b, dk, db = graphs.lsqm(np.log(freq[3:]), np.log(beta[3:]))
# ax.plot(np.log(freq[3:]), k * np.log(freq[3:]) + b)
# ax.errorbar(np.log(freq), np.log(beta), dbeta / beta, 0, '.')
# plt.show()

# print(f'k2 = {k:.3f} +- {dk:.3f}')
