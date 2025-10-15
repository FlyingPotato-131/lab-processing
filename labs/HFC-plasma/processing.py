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

height = 37.4 #mm

#plot clean
fig, ax = graphs.basePlot()

for image_path in sorted(os.listdir("clean_crop")):
	# Загружаем изображение
	img = Image.open(f"clean_crop/{image_path}").convert("RGB")
	# print(image_path)

	# Переводим в массив numpy
	img_array = np.array(img)

	r_channel_avg = img_array[:, :, 0].mean(axis=1)[::-1]
	g_channel_avg = img_array[:, :, 1].mean(axis=1)[::-1]
	b_channel_avg = img_array[:, :, 2].mean(axis=1)[::-1]

	# Усреднённая яркость (среднее от суммы)
	sum_channel_avg = (r_channel_avg + g_channel_avg + b_channel_avg) / 3

	x = np.arange(0, sum_channel_avg.size, 1) #x axis for plotting

	ax.plot(x / sum_channel_avg.size * height, sum_channel_avg / np.max(sum_channel_avg), linestyle='-', label=f"Среднее (R+G+B)/3, {image_path[:-4]} torr")

plt.legend()
plt.xlabel("z, pix")
plt.ylabel("intensity")
plt.title("intensity distribution")
plt.show()

#plot dirty
fig, ax = graphs.basePlot()

for image_path in sorted(os.listdir("dirty_crop")):
	# Загружаем изображение
	img = Image.open(f"dirty_crop/{image_path}").convert("RGB")
	# print(image_path)

	# Переводим в массив numpy
	img_array = np.array(img)

	r_channel_avg = img_array[:, :, 0].mean(axis=1)[::-1]
	g_channel_avg = img_array[:, :, 1].mean(axis=1)[::-1]
	b_channel_avg = img_array[:, :, 2].mean(axis=1)[::-1]


	# Усреднённая яркость (среднее от суммы)
	sum_channel_avg = (r_channel_avg + g_channel_avg + b_channel_avg) / 3
	ax.plot(sum_channel_avg / np.max(sum_channel_avg), linestyle='-', label=f"Среднее (R+G+B)/3, {image_path[:-4]} torr")

plt.legend()
plt.xlabel("z, pix")
plt.ylabel("intensity")
plt.title("intensity distribution")
plt.show()
