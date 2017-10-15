# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 20:27:41 2017

@author: naze
"""

import sys
import os
from time import sleep
import nibabel as nib
import numpy as np
from scipy.spatial.distance import cdist
import matplotlib.pylab as plt
from scipy.sparse import csr_matrix, lil_matrix
from scipy.io import mmwrite

#PRD = os.environ['PRD']
PRD = '/usr/local/freesurfer/subjects/R005484092'
mesh_size = 10242  # mesh size per hemisphere, to be given as argument
os.chdir(os.path.join(PRD, 'connectivity', 'ascii_tck'))
tck_file_list = os.listdir(os.path.join(PRD, 'connectivity', 'ascii_tck'))
tck_file_list = sorted(tck_file_list, key=str.lower)
print('extracting beginning and end coordinates of '+repr(len(tck_file_list))+' streamlines... and generate custom size connectivity matrix')

# get transform to get from dti to mri reference space
mri_img = nib.load(os.path.join(PRD, 'mri', 'T1.mgz'))
dti_img = nib.load(os.path.join(PRD, 'data', 'DWI', 'R005484092_visit_1_DTI_S152924.nii.gz'))
dti_vox2mri_vox = np.linalg.inv(mri_img.affine).dot(dti_img.affine)

#tck_bounds = open(os.path.join(PRD, 'connectivity','tckbounds.txt'), 'w')
lh_vtx = np.fromfile(os.path.join(PRD,'surface','lh_vertices_low.txt'),dtype=float, count=-1, sep=' ')
rh_vtx = np.fromfile(os.path.join(PRD,'surface','rh_vertices_low.txt'),dtype=float, count=-1, sep=' ')
lh_vtx = lh_vtx.reshape(mesh_size,3)
rh_vtx = rh_vtx.reshape(mesh_size,3)
lh_rh_vtx = np.array([lh_vtx, rh_vtx]).reshape(2*mesh_size,3)

custom_connectivity = csr_matrix((2*mesh_size, 2*mesh_size), dtype=np.int32)  # creates sparse matrix, mesh_size is for one hemisphere, whole brain connectome uses both hemispheres
tract_lengths = lil_matrix((2*mesh_size, 2*mesh_size))
dists = np.zeros((2,len(tck_file_list))) # first, last

#d2s = np.loadtxt(os.path.join(PRD,'connectivity', 'diffusion_2_struct.mat'))
#M = d2s[:3,:3]
#A = d2s[:3,3]
#M = [[1,-0,0],[-0,1,0],[-0,-0,1]] 
#A = [-126, -113.5634307861328, -54.10147857666016]

for fname in tck_file_list:
    # get first and last lines
    #with open(fname, "r") as f:    
    '''
    first = f.readline().rstrip()       # Read the first line.
    f.seek(-2, os.SEEK_END)             # Jump to the second last byte.
    while f.read(1) != b"\n":           # Until EOL is found...
        f.seek(-2, os.SEEK_CUR)         # ...jump back the read byte plus one more.
    last = f.readline().rstrip()        # Read last line.
    '''
    #lines = f.readlines()
    lines = np.loadtxt(fname)
    
    # transform into numpy arrays
    #first = np.fromstring(lines[0], dtype=np.float, sep=' ' )    
    #last = np.fromstring(lines[-1], dtype=np.float, sep=' ' )    
    #first = np.dot(M,lines[0].T)+A
    #last = np.dot(M,lines[-1].T)+A
    first = lines[0]
    last = lines[-1]
    
    # go from dti to mri space
    #first = nib.affines.apply_affine(dti_vox2mri_vox, first)
    #last = nib.affines.apply_affine(dti_vox2mri_vox, last)

    # find nearest neigbour vertices on cortical mesh
    first_dist = cdist(lh_rh_vtx, np.array([first]), 'euclidean')
    first_min_idx = np.argmin(first_dist)
    
    last_dist = cdist(lh_rh_vtx, np.array([last]), 'euclidean')
    last_min_idx = np.argmin(last_dist)

    # update connectivity matrix
    #if (first_min_idx != last_min_idx):
    custom_connectivity[last_min_idx, first_min_idx] += 1
    tract_lengths[last_min_idx, first_min_idx] += len(lines)
    dists[:,tck_file_list.index(fname)] = [np.min(first_dist), np.min(last_dist)]
    sys.stdout.write('\r')        
    sys.stdout.write(fname+' done')
    #print [np.min(first_dist), np.min(last_dist)]
    sys.stdout.flush()
    #sleep(0.25)

# post-processing connectivity tract lengths to get average tract length between 2 nodes
tract_lengths /= custom_connectivity
#custom_connectivity += custom_connectivity.T 
#custom_connectivity *= np.eye(mesh_size)/2

#np.savetxt('custom_connectivity.txt', custom_connectivity)
os.chdir(os.path.join(PRD, 'connectivity'))
mmwrite('sparse_custom_connectivity', custom_connectivity)  # saves as text file

c_fig = plt.figure()
plt.spy(custom_connectivity)
plt.show()

d_fig = plt.figure()
plt.subplot(1,2,1)
plt.hist(dists)
plt.subplot(1,2,2)
plt.hist(tract_lengths)
plt.show()

    
