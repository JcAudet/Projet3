import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from qutip import *

Data = np.loadtxt('Data.dat')
X = np.loadtxt('Domain_dim0.dat')
Y = np.loadtxt('Domain_dim1.dat')

X,Y = np.meshgrid(X,Y)

dat_max = np.abs(Data).max()

fig, ax = plt.subplots()

c = ax.pcolormesh(X, Y, Data)

fig.colorbar(c, ax=ax)

plt.show()

#fig = plt.figure()
#ax = fig.add_subplot(111,projection='3d')

#ax.plot_surface(X,Y,Data)

#plt.show()

