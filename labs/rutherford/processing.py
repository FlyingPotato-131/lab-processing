import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, pi

foil0 = csvreader.readData("0deg-foil.csv")
foil10 = csvreader.readData("10deg-foil.csv")
foil25 = csvreader.readData("25deg-foil.csv")
empty0 = csvreader.readData("0deg.csv")
empty10 = csvreader.readData("10deg.csv")
empty25 = csvreader.readData("25deg.csv")

midchan = np.sum(empty0[1750:1860, 0] * np.arange(1750, 1860, 1)) / np.sum(empty0[1750:1860, 0])
dmidchan = sqrt(np.sum((empty0[1750:1860, 0] * (np.arange(1750, 1860, 1) - midchan) / np.sum(empty0[1750:1860, 0]))**2))
print(f"центральный канал {midchan:.0f} +- {dmidchan:.0f}")
calib = 5.15 / midchan * 1e3
dcalib = 5.15 / midchan**2 * dmidchan * 1e3
print(f"калибровочный коэффициент ({calib:.3f} +- {dcalib:.3f}) кэв/кан")
print()

fig, ax = graphs.basePlot()
ax.plot(empty0, label="без фольги")
ax.plot(foil0, label="с фольгой")

print(f"n_0 = {np.sum(empty0[1750:1860, 0]):.0f} +- {sqrt(np.sum(empty0[1760:1860, 0])):.0f}")
print(f"n_foil0 = {np.sum(foil0[1750:1860, 0]):.0f} +- {sqrt(np.sum(foil0[1760:1860, 0])):.0f}")

midfoil0 = np.sum(foil0[1750:1860, 0] * np.arange(1750, 1860, 1)) / np.sum(foil0[1750:1860, 0])
dmidfoil0 = sqrt(np.sum((foil0[1750:1860, 0] * (np.arange(1750, 1860, 1) - midfoil0) / np.sum(foil0[1750:1860, 0]))**2))
print(f"Δε = ({(midchan - midfoil0) * calib:.1f} +- {(dmidchan + dmidfoil0) * calib + (midchan - midfoil0) * dcalib:.1f}) кэв")

print()

plt.xlabel("канал")
plt.ylabel("количество зарегестрированных частиц")
plt.title("спектр зарегестрированных частиц при φ=0deg")
plt.legend()
plt.show()

fig, ax = graphs.basePlot()
ax.plot(empty10, label="без фольги")
ax.plot(foil10, label="с фольгой")

print(f"n_10 = {np.sum(empty10[1750:1860, 0]):.0f} +- {sqrt(np.sum(empty10[1760:1860, 0])):.0f}")
print(f"n_foil10 = {np.sum(foil10[1750:1860, 0]):.0f} +- {sqrt(np.sum(foil10[1760:1860, 0])):.0f}")

midempty10 = np.sum(empty10[1750:1860, 0] * np.arange(1750, 1860, 1)) / np.sum(empty10[1750:1860, 0])
dmidempty10 = sqrt(np.sum((empty10[1750:1860, 0] * (np.arange(1750, 1860, 1) - midempty10) / np.sum(empty10[1750:1860, 0]))**2))
midfoil10 = np.sum(foil10[1750:1860, 0] * np.arange(1750, 1860, 1)) / np.sum(foil10[1750:1860, 0])
dmidfoil10 = sqrt(np.sum((foil10[1750:1860, 0] * (np.arange(1750, 1860, 1) - midfoil10) / np.sum(foil10[1750:1860, 0]))**2))
print(f"Δε = ({(midempty10 - midfoil10) * calib:.1f} +- {(dmidempty10 + dmidfoil10) * calib + (midempty10 - midfoil10) * dcalib:.1f}) кэв")

print()

plt.xlabel("канал")
plt.ylabel("количество зарегестрированных частиц")
plt.title("спектр зарегестрированных частиц при φ=10deg")
plt.legend()
plt.show()

fig, ax = graphs.basePlot()
ax.plot(empty25, label="без фольги")
ax.plot(foil25, label="с фольгой")

print(f"n_25 = {np.sum(empty25[1750:1860, 0]):.0f} +- {sqrt(np.sum(empty25[1760:1860, 0])):.0f}")
print(f"n_foil25 = {np.sum(foil25[1750:1860, 0]):.0f} +- {sqrt(np.sum(foil25[1760:1860, 0])):.0f}")

midempty25 = np.sum(empty25[1750:1860, 0] * np.arange(1750, 1860, 1)) / np.sum(empty25[1750:1860, 0])
dmidempty25 = sqrt(np.sum((empty25[1750:1860, 0] * (np.arange(1750, 1860, 1) - midempty25) / np.sum(empty25[1750:1860, 0]))**2))
midfoil25 = np.sum(foil25[1750:1860, 0] * np.arange(1750, 1860, 1)) / np.sum(foil25[1750:1860, 0])
dmidfoil25 = sqrt(np.sum((foil25[1750:1860, 0] * (np.arange(1750, 1860, 1) - midfoil25) / np.sum(foil25[1750:1860, 0]))**2))
print(f"Δε = ({(midempty25 - midfoil25) * calib:.1f} +- {(dmidempty25 + dmidfoil25) * calib + (midempty25 - midfoil25) * dcalib:.1f}) кэв")

print()

plt.xlabel("канал")
plt.ylabel("количество зарегестрированных частиц")
plt.title("спектр зарегестрированных частиц при φ=25deg")
plt.legend()
plt.show()
