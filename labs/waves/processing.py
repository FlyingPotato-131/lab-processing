import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import math
import numpy as np
import matplotlib.pyplot as plt

import os
import scipy
import csv

directory = os.listdir(".") #get filenames
for file in directory:
	if file.endswith(".csv"): #select all .csv
		data = np.empty(0)
		with open(file, newline = '') as raw:
			reader = csv.reader(raw) #read file
			data = np.array(list(reader)[2:])[:, :3].astype('f') #remove title and trailing ''

		ch1_raw = scipy.fft.fft(data[:, 1]) #fourier transform of data, e^(-2 pi i * (k n) / N)
		ch2_raw = scipy.fft.fft(data[:, 2])
		N = ch1_raw.size

		fig, ax = graphs.basePlot() #plot raw fourier
		ax.plot(ch1_raw, label = "channel 1")
		ax.plot(ch2_raw, label = "channel 2")
		plt.title("Raw Fourier transform at freq " + file[:-4] + " Hz")
		plt.xlabel("k, ν = k / (Nτ)")
		plt.ylabel("a")
		plt.legend()
		plt.show()

		cutoff = 20

		ch1 = ch1_raw[:cutoff] #remove frequencies with k >= 20 as it is noise
		ch2 = ch2_raw[:cutoff]

		print(f'initial signal freq f = {file[:-4]}')
		print(f'ch1 base freq f = {np.argmax(np.abs(ch1)) / N / 2e-3:.2f} Hz') #calculate base frequency of signals
		print(f'ch2 base freq f = {np.argmax(np.abs(ch2)) / N / 2e-3:.2f} Hz')

		fig, ax = graphs.basePlot() #plot truncated fourier
		ax.plot(ch1, label = "channel 1")
		ax.plot(ch2, label = "channel 2")
		plt.title("Truncated (k < 20) Fourier transform at freq " + file[:-4] + " Hz")
		plt.xlabel("k, ν = k / (Nτ)")
		plt.ylabel("a")
		plt.legend()
		plt.show()

		nu0 = 0.417 #t = n / N / nu0
		time = np.arange(0, N, 1) / N / nu0 #time axis

		def ch1_re(t): #function of restored signal
			return 1 / N * np.sum(ch1 * np.exp(2 * np.pi * 1j * np.arange(0, cutoff, 1) / N / 2e-3 * t))

		def ch2_re(t):
			return 1 / N * np.sum(ch2 * np.exp(2 * np.pi * 1j * np.arange(0, cutoff, 1) / N / 2e-3 * t))

		scale = 3

		fig, ax = graphs.basePlot() #plot restored signal
		ax.plot(time, data[:, 1], 'y', linewidth = 0.5, label = "channel 1")
		ax.plot(time, data[:, 2], 'tab:orange', linewidth = 0.5, label = "channel 2")
		ax.plot(time, [scale * ch1_re(t) for t in time], 'r', linewidth = 2, label = "channel 1 restored")
		ax.plot(time, [scale * ch2_re(t) for t in time], 'k', linewidth = 2, label = "channel 2 restored")
		plt.title("Signal and restored signal at freq " + file[:-4] + " Hz")
		plt.xlabel("t, s")
		plt.ylabel("V, V")
		plt.legend()
		plt.show()
