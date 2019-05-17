from qutip import *
import numpy as np

a = np.array([1,2,3,4,5,6,7,8,9])

print(a)

a.transpose()

print(a)

M = vec2mat(a)

print(M)

