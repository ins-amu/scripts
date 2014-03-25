import os
import numpy as np
from soma import aims
PRD = os.environ['PRD']
os.chdir(os.path.join(PRD, 'surface'))
mesh = aims.AimsTimeSurface( 3 )
# a mesh has a header
mesh.header()[ 'cortex_high' ] = 'cortex_high'
vert = mesh.vertex()
poly = mesh.polygon()
c = np.loadtxt('rh_vertices_high.txt')
vert.assign( [ aims.Point3df( x ) for x in c ] )
pol = np.loadtxt('rh_triangles_high.txt')
poly.assign( [ aims.AimsVector(x, dtype='U32',dim=3) for x in pol ] )
# write result
aims.write( mesh, 'rh_mesh_high.mesh' )