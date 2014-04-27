import os
import numpy as np
PRD = os.environ['PRD']
os.chdir(os.path.join(PRD, 'surface'))

with open('rh.pial.asc', 'r') as f:
    f.readline()
    nb_vert = f.readline().split(' ')[0]
    read_data = [[np.double(line.rstrip('\n').split()[0]),
                 np.double(line.rstrip('\n').split()[1]),
                 np.double(line.rstrip('\n').split()[2])] for line in f]

a = np.array(read_data)
vert_high = a[0:int(nb_vert), 0:3]
np.savetxt('rh_vertices_high.txt', vert_high)
tri_high = a[int(nb_vert):, 0:3]
np.savetxt('rh_triangles_high.txt', tri_high)
