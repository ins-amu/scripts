import sys
import numpy as np

vert = np.loadtxt(sys.argv[1])
tri = np.loadtxt(sys.argv[2])

with open(sys.argv[3], 'w') as f:
    f.write('OFF\n') 
    f.write('{} {} {}\n'.format(int(vert.shape[0]), int(tri.shape[0]), 0)) 

with open(sys.argv[3], 'ab') as f:
    np.savetxt(f, vert, fmt='%.6f')
    np.savetxt(f, np.hstack([np.ones((tri.shape[0],1))*3, tri]), fmt='%d')
    
