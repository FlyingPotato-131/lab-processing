from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
import scipy
import os

import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

# import seaborn as sns

def matrix(content): #code stolen from builtin plot.py
    # r - radial
    xRaw = content[2]
    xArrRaw = xRaw.split()
    xArr = []
    for i in range(len(xArrRaw)):
        xArr.append(float(xArrRaw[i]))
    x = np.array(xArr)

    # z and matrix
    z = []
    without_header = content[3:]
    data = []
    for line in without_header:
        data.append(line.strip())

    m = []
    for line in data:
        arr_raw = line.split()
        arr = arr_raw[1:]
        zOne = float(arr_raw[0])
        z.append(zOne)
        for i in range(len(arr)):
            arr[i] = float(arr[i])
        m.append(arr)
    return x, np.array(z), np.array(m)

index = ["N, 5torr", "N 2.5torr", "N 1torr", "Ar 5torr", "Ar 2.5torr", "Ar 1torr"] #index experiments

fig1, ax1 = graphs.basePlot()
fig2, ax2 = graphs.basePlot()

for image_path in sorted(os.listdir("images")):
    # Загружаем изображение
    img = Image.open(f"images/{image_path}").convert("RGB")
    # print(image_path)

    # Переводим в массив numpy
    img_array = np.array(img)

    r_channel_avg = img_array[:, :, 0].mean(axis=0)[::-1]
    g_channel_avg = img_array[:, :, 1].mean(axis=0)[::-1]
    b_channel_avg = img_array[:, :, 2].mean(axis=0)[::-1]


    # Усреднённая яркость (среднее от суммы)
    sum_channel_avg = (r_channel_avg + g_channel_avg + b_channel_avg) / 3

    z = np.arange(0, sum_channel_avg.shape[0], 1)
    def gauss(x, A, mu, sigma):
        return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

    params = scipy.optimize.curve_fit(gauss, z, sum_channel_avg) #подгоняем нормальное распределение
    height = 260
    diam = 25.3 #external tube diameter, mm

    # if(int(image_path[0]) == 1 or int(image_path[0]) == 4):
    #     # plt.figure(figsize=(12, 6)) #create graph
    #     fig, ax = graphs.basePlot()

    #open theoretical data (if exists)
    # if os.path.isdir(f"theor/{image_path[0]}"):
    #     with open(f"theor/{image_path[0]}/Temp.txt", encoding = "cp1251") as dist_file: #read file
    #         dist_theor = dist_file.read().splitlines()

    #     r, z_t, data = matrix(dist_theor)        
    #     avgdata = data.mean(axis=1)
    #     # print(data)

    #     plt.plot(z_t * 10, avgdata / np.max(avgdata) * np.max(sum_channel_avg)) #plot converting z to mm and normalizing units

    # if os.path.isdir(f"theor/{image_path[0]}"):
        # print(image_path[0])
        # dist_theor = csvreader.readTable(f"theor/{image_path[0]}/Power.txt", 251, titleSize = 3, split = ' ', enc = "cp1251")

    # Строим
    if(int(image_path[0]) < 4):
        ax1.plot(z / height * diam, sum_channel_avg / np.max(sum_channel_avg), linestyle='-', label=f"Среднее (R+G+B)/3, {index[int(image_path[0]) - 1]}")
        ax1.plot(z / height * diam, gauss(z, *(params[0])) / np.max(sum_channel_avg), linestyle = '--', label = f"подгон нормального распределения, {index[int(image_path[0]) - 1]}")
        ax1.plot([params[0][1] / height * diam, params[0][1] / height * diam], [0, 1], linestyle = "--", label = f"max излучения z = {params[0][1] / height * diam:.0f} mm, {index[int(image_path[0]) - 1]}")
    else:    
        ax2.plot(z / height * diam, sum_channel_avg / np.max(sum_channel_avg), linestyle='-', label=f"Среднее (R+G+B)/3, {index[int(image_path[0]) - 1]}")
        ax2.plot(z / height * diam, gauss(z, *(params[0])) / np.max(sum_channel_avg), linestyle = '--', label = f"подгон нормального распределения, {index[int(image_path[0]) - 1]}")
        ax2.plot([params[0][1] / height * diam, params[0][1] / height * diam], [0, 1], linestyle = "--", label = f"max излучения z = {params[0][1] / height * diam:.0f} mm, {index[int(image_path[0]) - 1]}")
    # plt.title(f"Усреднённые интенсивности каналов RGB, {index[int(image_path[0]) - 1]}")
plt.xlabel("z, mm")
plt.ylabel("Интенсивность")
ax1.legend()
ax2.legend()
    # plt.show()
    # if(int(image_path[0]) == 3 or int(image_path[0]) == 6):
plt.show()