import sys

sys.path.append("../../")
import csvreader
import graphs

data1 = csvreader.readTable('data-e.csv', 2, titleSize = 0)
data2 = csvreader.readTable('data2.csv', 2, titleSize = 0)
data3 = csvreader.readTable('data3.csv', 2, titleSize = 0)
data4 = csvreader.readTable('data4.csv', 2, titleSize = 1)

graphs.plot(data1[0:, 1], data1[0:, 0], title = '346 - 322', ylabel = 'F, n', xlabel = 'dx, 0.1*mm')
graphs.plot(data2[0:, 1], data2[0:, 0], title = '177 - 490', ylabel = 'F, n', xlabel = 'dx, 0.1*mm')
graphs.plot(data3[0:, 1], data3[0:, 0], title = '240 - 430', ylabel = 'F, n', xlabel = 'dx, 0.1*mm')
graphs.plot(data4[0:, 1], data4[0:, 0], title = '492 - 172', ylabel = 'F, n', xlabel = 'dx, 0.1*mm')

#погрешность есть половина цены деления
