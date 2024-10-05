import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, pi

data = csvreader.readData("data.csv")

calib = np.array([[1.17, data[0,0]],[1.33, data[2,0]],[0.662, data[2,1]],[1.274, data[2,2]]])

K, B, dK, dB = graphs.plotlsqm(calib[:,0], calib[:,1], title = "калибровочный график", xlabel = "энергия фотона, МэВ", ylabel = "номер канала")

print(K, B)

data[0,:] = (data[0,:] - B) / K
data[1,:] = (data[1,:]) / K
data[2,:] = (data[2,:] - B) / K
data[3,:] = (data[3,:]) / K
data[4,:] = (data[4,:] - B) / K

print("приведенная к МэВ таблица результатов")
print(data)
print()

Ri = data[3,:] / data[2,:]
dRi = np.sqrt((0.1 / K / data[2,:])**2 + (data[3,:] / data[2,:] * 0.1 / K))
print("энергетическое разрешение")
print(Ri)
print("+-", dRi)
print()

mc2 = 0.511
Etheor = np.array(data[2,:] / (1 + mc2 / 2 / data[2,:]))
graphs.plot(Etheor, data[4,:], title = "экспериментальные и рассчетные границы Комптоновского излучения", xlabel = "рассчетная энергия, МэВ", ylabel = "экспериментальная энергия, МэВ")

graphs.plotlsqm(1 / np.concatenate([data[2,:3], [data[2, 4]]]), np.concatenate([Ri[:3], [Ri[4]]])**2, 0.1 / K / np.concatenate([data[2,:3], [data[2, 4]]])**2, 2 * np.concatenate([Ri[:3], [Ri[4]]]) * np.concatenate([dRi[:3], [dRi[4]]]), title = "проверка зависимости (6) (кроме Am)", xlabel = "1 / E, MeV^-1", ylabel = "Ri^2")

cs137 = csvreader.readData("Cs137.csv")
cs137alt = csvreader.readData("Cs137-alt.csv")

fig, ax = graphs.basePlot()
ax.plot(cs137[:, 0], cs137[:, 1], label = "our")
ax.plot(cs137alt[:, 0], cs137alt[:, 1], label = "alt")
plt.legend()
plt.title("Спектры цезия 137 на нашей и соседней установках")
plt.xlabel("канал")
plt.ylabel("количество фотонов")
plt.show()

Am241 = csvreader.readData("Am241.csv")
fig, ax = graphs.basePlot()
ax.plot(Am241[:, 0], Am241[:, 1])
plt.title("Спектры Америция 241")
plt.xlabel("канал")
plt.ylabel("количество фотонов")
plt.show()

Co60 = csvreader.readData("Co60.csv")
fig, ax = graphs.basePlot()
ax.plot(Co60[:, 0], Co60[:, 1])
plt.title("Спектры Кобальта 60")
plt.xlabel("канал")
plt.ylabel("количество фотонов")
plt.show()

Eu152 = csvreader.readData("Eu152.csv")
fig, ax = graphs.basePlot()
ax.plot(Eu152[:, 0], Eu152[:, 1])
plt.title("Спектры Европия 152")
plt.xlabel("канал")
plt.ylabel("количество фотонов")
plt.show()

Na22 = csvreader.readData("Na22.csv")
fig, ax = graphs.basePlot()
ax.plot(Na22[:, 0], Na22[:, 1])
plt.title("Спектры Натрия 22")
plt.xlabel("канал")
plt.ylabel("количество фотонов")
plt.show()