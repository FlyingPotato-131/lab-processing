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

index = ["", "", "Ar 2.5torr", "Ar 1torr", "", "Ar 5torr"]

for image_path in os.listdir("images"):
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

    #open theoretical data
    # with open(f"theor/{image_path[0]}/Power.txt") as 
    if os.path.isdir(f"theor/{image_path[0]}"):
        print(image_path[0])
        dist_theor = csvreader.readTable(f"theor/{image_path[0]}/Power.txt", 251, titleSize = 3, split = ' ', enc = "cp1251")

    # Строим
    plt.figure(figsize=(12, 6))
    plt.plot(z / height * diam, sum_channel_avg, color='black', linestyle='-', label="Среднее (R+G+B)/3")
    plt.plot(z / height * diam, gauss(z, *(params[0])), color = 'red', linestyle = '--', label = "подгон нормального распределения")
    plt.plot([params[0][1] / height * diam, params[0][1] / height * diam], [0, np.max(sum_channel_avg)], color = "violet", linestyle = "--")

    plt.title(f"Усреднённые интенсивности каналов RGB, {index[int(image_path[0]) - 1]}")
    plt.xlabel("z, mm")
    plt.ylabel("Интенсивность")
    plt.legend()
    plt.grid(True)
    plt.show()