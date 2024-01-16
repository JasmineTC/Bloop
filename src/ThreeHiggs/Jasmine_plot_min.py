import numpy as np
import matplotlib.pylab as plt

with open('results_test.txt') as f:
    data = f.readlines()

data_array = np.zeros(len(data))
for i in range(len(data)):
    data_point = data[i][0] + data[i][1] + data[i][2] + data[i][3]
    data_array[i] = float(data_point)

plt.plot(data_array,'.')
plt.title('Global min')
plt.xlabel('index')
plt.ylabel('vev')
plt.show()


