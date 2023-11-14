import csvreader
import graphs
import numpy as np

data = csvreader.readTable("labs/spectre-analysis/data/sq-100us/3000hz.csv", 2, titleSize = 3)
np.random.shuffle(data)
graphs.plot(data[0:, 0], data[0:, 1], plotFmt = '-')

