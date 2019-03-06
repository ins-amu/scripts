import numpy as np
import os
import sys

name_file = sys.argv[1]
PRD = os.environ['PRD']
FS = os.environ['FS']
SUBJ_ID = os.environ['SUBJ_ID']
os.chdir(os.path.join(PRD, 'surface'))

with open(os.path.join(FS, SUBJ_ID, 'bem', name_file + '.asc'), 'r') as f:
    f.readline()
    nb_vert = f.readline().split(' ')[0]
    read_data = [[np.double(line.rstrip('\n').split()[0]),
                 np.double(line.rstrip('\n').split()[1]),
                 np.double(line.rstrip('\n').split()[2])] for line in f]

a = np.array(read_data)
vert = a[0:int(nb_vert), 0:3]
np.savetxt(os.path.join(PRD, SUBJ_ID, 'surface', name_file +'_vertices.txt'), vert, fmt='%.6f %.6f %.6f')
tri = a[int(nb_vert):, 0:3]
np.savetxt(os.path.join(PRD, SUBJ_ID, 'surface', name_file +'_triangles.txt'), tri, fmt='%d %d %d')
