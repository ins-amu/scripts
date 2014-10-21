import os
import sys
import numpy as np
from soma import aims
PRD = os.environ['PRD']
os.chdir(os.path.join(PRD, 'surface'))
rl = sys.argv[1]

mesh = aims.AimsTimeSurface( 3 )
# a mesh has a header
mesh.header()[ 'cortex_high' ] = 'cortex_high'
vert = mesh.vertex()
poly = mesh.polygon()
c = np.loadtxt(rl+'_vertices_high.txt')
vert.assign( [ aims.Point3df( x ) for x in c ] )
pol = np.loadtxt(rl+'_triangles_high.txt')
poly.assign( [ aims.AimsVector(x, dtype='U32',dim=3) for x in pol ] )
# write result
aims.write( mesh, rl+'_mesh_high.mesh' )
