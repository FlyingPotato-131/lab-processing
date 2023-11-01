import sys

sys.path.append("../../")
import csvreader
import graphs
import calculations

data = csvreader.readData('experiment-data.csv')
print(data)