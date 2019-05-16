import matplotlib.pyplot as plt
import numpy as np

Data = np.loadtxt('test.dat')

for x in Data:
    print(x)

plt.plot(Data,Data,label='linear')
plt.plot(Data,Data**2,label='Quadratic')

plt.legend()


plt.show()




