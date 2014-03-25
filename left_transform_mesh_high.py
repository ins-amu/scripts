import os
from soma import aims
import numpy
PRD = os.environ['PRD']
os.chdir(os.path.join(PRD, 'surface'))
mesh = aims.AimsTimeSurface( 3 )
# a mesh has a header
mesh.header()[ 'cortex_high' ] = 'cortex_high'
vert = mesh.vertex()
poly = mesh.polygon()
c = numpy.loadtxt('lh_vertices_high.txt')
vert.assign( [ aims.Point3df( x ) for x in c ] )
pol = numpy.loadtxt('lh_triangles_high.txt')
poly.assign( [ aims.AimsVector(x, dtype='U32',dim=3) for x in pol ] )
# write result
aims.write( mesh, 'lh_mesh_high.mesh' )
