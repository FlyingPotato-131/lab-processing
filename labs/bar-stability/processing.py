import sys

sys.path.append("../../")
import csvreader
import graphs

data = csvreader.readTable('data3.csv', 2, titleSize = 0)

graphs.plot(data[0:, 1], data[0:, 0])

#погрешность есть половина цены деления