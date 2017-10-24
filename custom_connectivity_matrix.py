# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 20:27:41 2017

@author: naze
"""

import sys
import os
import nibabel as nib
import numpy as np
from scipy.spatial.distance import cdist
import matplotlib.pylab as plt
from scipy.sparse import csr_matrix, lil_matrix
from scipy.io import mmwrite

#PRD = os.environ['PRD']
PRD = '/usr/local/freesurfer/subjects/R005484092_2'
mesh_size = 10242  # mesh size per hemisphere, to be given as argument
#os.chdir(os.path.join(PRD, 'connectivity', 'tmp_ascii_tck'))
tck_file_list = os.listdir(os.path.join(PRD, 'connectivity', 'tmp_ascii_tck'))
tck_file_list = np.asarray(sorted(tck_file_list, key=str.lower))
print('extracting beginning and end coordinates of '+repr(len(tck_file_list))+' streamlines... and generate custom size connectivity matrix')

# get transforms to from diffusion to structural reference spaces
img = nib.load(os.path.join(PRD,'connectivity','brain.nii.gz'))
trans_vox2mm = img.affine
M_vox2mm = trans_vox2mm[:3, :3]
A_vox2mm = trans_vox2mm[:3,3]

brain = img.get_data()
M_RAS2LAS = [[-1,0,0],[0,1,0],[0,0,1]]
A_RAS2LAS = [brain.shape[0], 0, 0]

# get decimated cortical surface (LAS, mm)
lh_vtx = np.loadtxt(os.path.join(PRD,'surface','lh_vertices_low.txt'))
rh_vtx = np.loadtxt(os.path.join(PRD,'surface','rh_vertices_low.txt'))
vtx = np.vstack([lh_vtx, rh_vtx])

#custom_connectivity = csr_matrix((2*mesh_size, 2*mesh_size), dtype=np.int32)  # creates sparse matrix, mesh_size is for one hemisphere, whole brain connectome uses both hemispheres
#tract_lengths = lil_matrix((2*mesh_size, 2*mesh_size))
custom_connectivity = np.zeros((2*mesh_size, 2*mesh_size))  # creates sparse matrix, mesh_size is for one hemisphere, whole brain connectome uses both hemispheres

max_dist = 5 # max distance from the fiber tract to the cortical surface to be considered 
nsamples=100000
tract_lengths = []  #lil_matrix((2*mesh_size, 2*mesh_size))
dists = np.zeros((2,nsamples)) # first, last

#d2s = np.loadtxt(os.path.join(PRD,'connectivity', 'diffusion_2_struct.mat'))
#M = d2s[:3,:3]
#A = d2s[:3,3]
#M = [[1,-0,0],[-0,1,0],[-0,-0,1]] 
#A = [-126, -113.5634307861328, -54.10147857666016]
randidx = np.random.randint(len(tck_file_list), size=(1,nsamples))
randidx = np.sort(randidx).flatten()
for i in np.arange(0,nsamples):
    fname = tck_file_list[randidx[i]]
    '''
    first = f.readline().rstrip()       # Read the first line.
    f.seek(-2, os.SEEK_END)             # Jump to the second last byte.
    while f.read(1) != b"\n":           # Until EOL is found...
        f.seek(-2, os.SEEK_CUR)         # ...jump back the read byte plus one more.
    last = f.readline().rstrip()        # Read last line.
    '''
    #lines = f.readlines()
    lines = np.loadtxt(os.path.join(PRD,'connectivity','tmp_ascii_tck',fname))
    lines = np.dot(lines, M_RAS2LAS) + A_RAS2LAS
    lines = np.dot(lines, M_vox2mm) + A_vox2mm
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
    dists[:,i] = [np.min(first_dist), np.min(last_dist)]
    if np.max(dists[:,i]) < max_dist :
        custom_connectivity[last_min_idx, first_min_idx] += 1
        tract_lengths.append(len(lines))
    
    sys.stdout.write('\r')        
    sys.stdout.write(fname+' done')
    #print [np.min(first_dist), np.min(last_dist)]
    sys.stdout.flush()

# post-processing connectivity tract lengths to get average tract length between 2 nodes
#tract_lengths /= custom_connectivity
#custom_connectivity += custom_connectivity.T 
#custom_connectivity *= np.eye(mesh_size)/2

#np.savetxt('custom_connectivity.txt', custom_connectivity)
mmwrite(os.path.join(PRD,'connectivity','sparse_custom_connectivity'), csr_matrix(custom_connectivity))  # saves as text file
#mmwrite('sparse_tract_lengths', csr_matrix(tract_lengths))  # saves as text file
#np.save('distances_between dti_and_mri_vertices', dists)  # saves as text file

c_fig = plt.figure()
plt.spy(csr_matrix(custom_connectivity), marker='.', markersize=0.5)
plt.show()

d_fig = plt.figure()
plt.subplot(1,2,1)
plt.hist(dists.flatten(), bins=np.arange(0,60,0.25))
plt.subplot(1,2,2)
plt.hist(tract_lengths, bins=np.arange(0,250,5))
plt.show()
 
