import os
import numpy as np
PRD = os.environ['PRD']
os.chdir(os.path.join(PRD, 'surface'))
a = np.loadtxt('rh_pial.txt') 
nb_vert =  np.loadtxt('rh_number_vertices_high.txt')
vert_high = a[0:int(nb_vert), 0:3]
np.savetxt('rh_vertices_high.txt', vert_high) 
tri_high = a[int(nb_vert):,0:3]
np.savetxt('rh_triangles_high.txt', tri_high) 