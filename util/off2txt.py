import numpy as np
import sys

name_file = sys.argv[1]

with open(name_file) as f:
    f.readline()
    num = f.readline().split(' ')

vert = np.loadtxt(name_file, skiprows=2, usecols=(0,1,2))  
vert = vert[:int(num[0]), :]
tri = np.loadtxt(name_file, skiprows=int(num[0])+2, usecols=(1,2,3))  

np.savetxt(sys.argv[2], vert, fmt='%.4f')
np.savetxt(sys.argv[3], tri, fmt='%d')

