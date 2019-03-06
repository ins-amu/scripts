import numpy as np
import os
import sys

rl = sys.argv[1]
PRD = os.environ['PRD']
os.chdir(os.path.join(PRD, 'surface'))

# read lh_info
with open(rl + 'info.txt', 'r') as f:
    lines = f.readlines()

c_ras_line = lines[32]
ista = c_ras_line.index(' (') + 2
iend = c_ras_line.index(')\n')
lc_ras = c_ras_line[ista:iend].split(',')
c_ras = np.array(lc_ras).astype('float')

with open(rl + '.pial.asc', 'r') as f:
    f.readline()
    nb_vert = f.readline().split(' ')[0]
    read_data = [[np.double(line.rstrip('\n').split()[0]),
                 np.double(line.rstrip('\n').split()[1]),
                 np.double(line.rstrip('\n').split()[2])] for line in f]

a = np.array(read_data)
vert_high = a[0:int(nb_vert), 0:3] + c_ras
np.savetxt(rl + '_vertices_high.txt', vert_high, fmt='%.6f %.6f %.6f')
tri_high = a[int(nb_vert):, 0:3]
np.savetxt(rl +'_triangles_high.txt', tri_high, fmt='%d %d %d')
