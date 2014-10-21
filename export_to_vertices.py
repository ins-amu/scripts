import os
import numpy as np
from soma import aims
import sys
rl = sys.argv[1]

PRD = os.environ['PRD']
os.chdir(os.path.join(PRD, 'surface'))
a = aims.read(rl + '_mesh_low.mesh')
b = np.zeros((len(a.vertex().list()),3))
g = a.vertex().list()
for i in range(len(g)):
    b[i,:] = g[i][0:3] 
np.savetxt(rl + '_vertices_low.txt', b)

a = aims.read(rl + '_mesh_low.mesh')
b = np.zeros((len(a.polygon().list()),3))
g = a.polygon().list()
for i in range(len(g)):
    b[i,:] = g[i][0:3] 
np.savetxt(rl + '_triangles_low.txt', b)
