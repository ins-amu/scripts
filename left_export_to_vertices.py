import os
import numpy as np
from soma import aims
PRD = os.environ['PRD']
os.chdir(os.path.join(PRD, 'surface'))
a = aims.read('lh_mesh_low.mesh')
b = np.zeros((len(a.vertex().list()),3))
g = a.vertex().list()
for i in range(len(g)):
    b[i,:] = g[i][0:3] 
np.savetxt('lh_vertices_low.txt', b)

a = aims.read('lh_mesh_low.mesh')
b = np.zeros((len(a.polygon().list()),3))
g = a.polygon().list()
for i in range(len(g)):
    b[i,:] = g[i][0:3] 
np.savetxt('lh_triangles_low.txt', b)