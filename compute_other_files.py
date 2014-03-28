import numpy as np
import os
PRD = os.environ['PRD']
SUBJ_ID = os.environ['SUBJ_ID']
from tvb.simulator.lab import *
default_cortex = surfaces.Cortex()
default_cortex.configure() 
orientations = default_cortex.region_orientation  
areas = default_cortex.region_areas  
centers = default_cortex.region_center 

# compute subcortical centers
corr_table = np.loadtxt('correspondance_table.txt')
for val in []:
	verts = np.load()
	curr_center = np.mean(verts, axis=0)
	indx = corr_table[np.nonzero(corr_table[:,0]==np.int(val)),2] - 1	
	centers[indx,:] = curr_center

np.savetxt(os.path.join(PRD, SUBJ_ID, 'connectivity/area'), areas, fmt='%.2f')
np.savetxt(os.path.join(PRD, SUBJ_ID, 'connectivity/orientation'), orientations, fmt='%.2f %.2f %.2f')

f = open('name_regions.txt', 'rb')
list_name = []
for line in f:
    list_name.append(line)
f.close()

f = open(os.path.join(PRD, SUBJ_ID, 'connectivity/position'), 'w')
for i, name in enumerate(list_name):
    f.write(str(name[:-1])+' ')
    for j in range(3):
        f.write('{:.4f} '.format(centers[i, j]))
    f.write('\n')
f.close()