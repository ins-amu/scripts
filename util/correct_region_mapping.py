import os
import sys
from copy import deepcopy
from collections import Counter
import numpy as np

rl = sys.argv[1]
PRD = os.environ['PRD']
region_mapping_corr = float(os.environ['REGION_MAPPING_CORR'])
os.chdir(os.path.join(PRD, 'surface'))


texture = np.loadtxt(rl + '_region_mapping_low_not_corrected.txt')
vert = np.loadtxt(rl + '_vertices_low.txt')
trian = np.loadtxt(rl + '_triangles_low.txt')
for _ in range(10):
    new_texture = deepcopy(texture)
    labels = np.unique(texture)
    for ilab in labels:
        iverts = (np.nonzero(texture==ilab)[0]).tolist()
        if len(iverts)>0:
            for inode in iverts:
                iall = trian[np.nonzero(trian==inode)[0]].flatten().tolist()
                #import pdb; pdb.set_trace()
                ineig = np.unique(list(filter(lambda x: x!=inode, iall))).astype('int')
                ivals = np.array(Counter(texture[ineig]).most_common()).astype('int')                
                if ivals[np.nonzero(ivals[:,0]==ilab), 1].shape[1]==0:
                    new_texture[inode] = Counter(texture[ineig]).most_common(1)[0][0]
                elif ivals[np.nonzero(ivals[:,0] == ilab), 1][0,0] < region_mapping_corr * len(ineig):
                    new_texture[inode] = Counter(texture[ineig]).most_common(1)[0][0]
    texture = deepcopy(new_texture) 

np.savetxt(rl + '_region_mapping_low.txt', new_texture)
