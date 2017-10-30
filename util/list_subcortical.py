import numpy as np
import os

PRD = os.environ['PRD']
os.chdir(os.path.join(PRD, 'surface', 'subcortical'))

for val in ['16','08','10','11','12','13','17','18','26','47','49','50','51','52','53','54','58']:
    val
    f = open('aseg_0'+val+'.srf', 'r')
    f.readline()
    data = f.readline()
    g = data.split(' ')
    nb_vert = int(g[0])
    nb_tri = int(g[1].split('\n')[0])
    f.close()
    a = np.loadtxt('aseg_0'+str(val)+'.srf', skiprows=2, usecols=(0,1,2))
    vert = a[:nb_vert]
    tri = a[nb_vert:].astype('int')
    np.savetxt('aseg_0'+str(val)+'_vert.txt', vert)
    np.savetxt('aseg_0'+str(val)+'_tri.txt', tri)
