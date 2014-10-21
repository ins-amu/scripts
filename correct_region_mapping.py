import os
import sys
rl = sys.argv[1]
PRD = os.environ['PRD']
region_mapping_corr = float(os.environ['region_mapping_corr'])
os.chdir(os.path.join(PRD, 'surface'))
from copy import deepcopy
from pylab import *
from collections import Counter

texture = loadtxt(rl + '_region_mapping_low_not_corrected.txt')
vert = loadtxt(rl + '_vertices_low.txt')
trian = loadtxt(rl + '_triangles_low.txt')
for _ in range(10):
    new_texture = deepcopy(texture)
    labels = np.unique(texture)
    #import pdb; pdb.set_trace()
    for ilab in labels:
        iverts = (np.nonzero(texture==ilab)[0]).tolist()
        if len(iverts)>0:
            for inode in iverts:
                iall = trian[np.nonzero(trian==inode)[0]].flatten().tolist()
                ineig = unique(filter(lambda x: x!=inode, iall)).astype('int')
                # import pdb; pdb.set_trace()
                ivals = np.array(Counter(texture[ineig]).most_common()).astype('int')                
                if ivals[np.nonzero(ivals[:,0]==ilab), 1].shape[1]==0:
                    new_texture[inode] = Counter(texture[ineig]).most_common(1)[0][0]
                elif ivals[np.nonzero(ivals[:,0] == ilab), 1][0,0] < region_mapping_corr * len(ineig):
                    new_texture[inode] = Counter(texture[ineig]).most_common(1)[0][0]
    texture = deepcopy(new_texture) 

savetxt(rl + '_region_mapping_low.txt', new_texture)
