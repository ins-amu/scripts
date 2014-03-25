import numpy as np
import os
PRD = os.environ['PRD']
os.chdir(os.path.join(PRD, 'surface'))
a = np.loadtxt('lh_pial.txt') 
nb_vert =  np.loadtxt('lh_number_vertices_high.txt')
vert_high = a[0:int(nb_vert), 0:3]
np.savetxt('lh_vertices_high.txt', vert_high) 
tri_high = a[int(nb_vert):,0:3]
np.savetxt('lh_triangles_high.txt', tri_high) 