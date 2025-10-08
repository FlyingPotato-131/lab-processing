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

index = ["N, 1torr", "Ar 2.5torr", "N 2.5torr", "N 5torr", "Ar 5torr", "Ar 5torr"] #index experiments

for image_path in os.listdir("images"):
    # Загружаем изображение
    img = Image.open(f"images/{image_path}").convert("RGB")
    print(image_path)

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

    plt.figure(figsize=(12, 6)) #create graph

    #open theoretical data (if exists)
    if os.path.isdir(f"theor/{image_path[0]}"):
        with open(f"theor/{image_path[0]}/Temp.txt", encoding = "cp1251") as dist_file: #read file
            dist_theor = dist_file.read().splitlines()

        r, z_t, data = matrix(dist_theor)        
        avgdata = data.mean(axis=1)
        # print(data)

        # plt.plot(z_t * 10, avgdata / np.max(avgdata) * np.max(sum_channel_avg)) #plot converting z to mm and normalizing units

    # if os.path.isdir(f"theor/{image_path[0]}"):
        # print(image_path[0])
        # dist_theor = csvreader.readTable(f"theor/{image_path[0]}/Power.txt", 251, titleSize = 3, split = ' ', enc = "cp1251")

    # Строим
    plt.plot(z / height * diam, sum_channel_avg, color='black', linestyle='-', label="Среднее (R+G+B)/3")
    plt.plot(z / height * diam, gauss(z, *(params[0])), color = 'red', linestyle = '--', label = "подгон нормального распределения")
    plt.plot([params[0][1] / height * diam, params[0][1] / height * diam], [0, np.max(sum_channel_avg)], color = "violet", linestyle = "--", label = f"max излучения z = {params[0][1] / height * diam:.0f} mm")

    plt.title(f"Усреднённые интенсивности каналов RGB, {index[int(image_path[0]) - 1]}")
    plt.xlabel("z, mm")
    plt.ylabel("Интенсивность")
    plt.legend()
    plt.grid(True)
    plt.show()