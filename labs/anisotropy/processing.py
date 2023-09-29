import sys

sys.path.append("../../")
import csvreader
import graphs

data = csvreader.readTable('data.csv', 4, titleSize = 1)

graphs.plot(data[0:, 0], data[0:, 2])
