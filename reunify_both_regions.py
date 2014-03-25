import os
PRD = os.environ['PRD']
SUBJ_ID = os.environ['SUBJ_ID']
os.chdir(os.path.join(PRD, 'surface'))
from pylab import *
lh_reg_map = np.loadtxt('lh_region_mapping_low.txt')
lh_vert = np.loadtxt('lh_vertices_low.txt')
lh_trian = np.loadtxt('lh_triangles_low.txt')
rh_reg_map = np.loadtxt('rh_region_mapping_low.txt')
rh_vert = np.loadtxt('rh_vertices_low.txt')
rh_trian = np.loadtxt('rh_triangles_low.txt')
vertices = vstack([lh_vert, rh_vert])
triangles = vstack([lh_trian,  rh_trian + lh_vert.shape[0]]) 
region_mapping = hstack([lh_reg_map, rh_reg_map])
np.savetxt(PRD+'/'+SUBJ_ID+'/surfaces/region_mapping.txt', region_mapping, fmt='%d', newline=" ")
np.savetxt(PRD+'/'+SUBJ_ID+'/surfaces/vertices.txt', vertices, fmt='%.2f')
np.savetxt(PRD+'/'+SUBJ_ID+'/surfaces/triangles.txt', triangles, fmt='%d %d %d')