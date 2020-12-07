import nibabel as nib
import os
import sys
import collections
import pdb
import numpy as np
from copy import deepcopy
from mpl_toolkits.mplot3d import Axes3D
from collections import Counter
import matplotlib.pyplot as plt


def surface_projection(SUBJ_ID, PRD, img, img_aff, label, subcortical, DISPLAY, 
	                   CHECK):
    """
    img: nifti of the parcellation used for projection
    img_aff: nfiti of the original parcellation from which to extract the affine
    subcortical: list of subcortical regions which are in the nifti
    """
    print('project original region mapping to the new parcellation')

    parcellation_new  = img.get_data()
    affine = img_aff.get_affine()
    verts_origin= np.loadtxt(os.path.join(PRD, SUBJ_ID, 'surface', 'vertices.txt')) 

    # project vertice coordinaates in voxel coordinaates in diffusion space
    verts_vox = np.round(np.dot(np.linalg.inv(affine), np.hstack([verts_origin, np.ones((verts_origin.shape[0], 1))]).T)[:3].T).astype(int)

    if CHECK=="yes" and len(DISPLAY)>0:
        # check that the surface and the parcellation are in the same space
        n = 100 
        slice_n = np.nonzero(verts_vox[:, 0]==n)
        plt.figure(figsize=(13, 5))
        plt.subplot(121)
        plt.scatter(verts_vox[slice_n, 1], verts_vox[slice_n, 2])
        plt.xlim((0, 250))
        plt.ylim((0, 250))
        plt.subplot(122)
        plt.imshow(parcellation_new[n, :, :].T, origin='lower', interpolation='nearest')
        plt.show()

    region_mapping = np.loadtxt(os.path.join(PRD, SUBJ_ID, 'region_mapping.txt')) 

    region_mapping_new = []
    for ivert, vert in enumerate(verts_vox):
        if region_mapping[ivert]==0: # region 0: keep the same
            region_mapping_new.append(0)
        else: # search for the closest region with a label
            # if ivert%100==0:print(ivert)
            vox_list, vox_list_add = [], []
            vox_list.append((vert[0], vert[1], vert[2]))
            vert_assign = [0]
            icount = 0
            while (np.unique(vert_assign)==np.array([0])).all(): # while a region mapping is not found
                icount += 1
                vox_list.extend(vox_list_add)
                for ivox in vox_list: # find the region mapping value for the voxes list
                    vert_assign.append(parcellation_new[ivox[0], ivox[1], ivox[2]])
                for ivox in vox_list: # add the nearest neighbout vox
                    if (ivox[0]+1, ivox[1], ivox[2]) not in vox_list + vox_list_add:
                        vox_list_add.append((ivox[0]+1, ivox[1], ivox[2]))
                    if (ivox[0]-1, ivox[1], ivox[2]) not in vox_list + vox_list_add:
                        vox_list_add.append((ivox[0]-1, ivox[1], ivox[2]))
                    if (ivox[0], ivox[1]+1, ivox[2]) not in vox_list + vox_list_add:
                        vox_list_add.append((ivox[0], ivox[1]+1, ivox[2]))
                    if (ivox[0], ivox[1]-1, ivox[2]) not in vox_list + vox_list_add:
                        vox_list_add.append((ivox[0], ivox[1]-1, ivox[2]))
                    if (ivox[0], ivox[1], ivox[2]+1) not in vox_list + vox_list_add:
                        vox_list_add.append((ivox[0], ivox[1], ivox[2]+1))
                    if (ivox[0], ivox[1], ivox[2]-1) not in vox_list + vox_list_add:
                        vox_list_add.append((ivox[0], ivox[1], ivox[2]-1))
                if icount>10:
                    raise ValueError('had too look more than 10 voxels away, \
                                      check parcellation/surface coregistration')
            final_vox_list = [i for i in vert_assign if i!=0]
            most_common_vox = collections.Counter(final_vox_list).most_common()[0][0]-1
            if len(subcortical)>0: # assign the mapping of subcortical regions to 0
                                   # as the region mapping only applies to the surface
                if most_common_vox in subcortical: most_common_vox=0
            region_mapping_new.append(most_common_vox)
    # add the subcortical regions at the end of the region mapping (TVB convention)
    # region_mapping_new.extend(subcortical) # depreciated?

    np.savetxt(os.path.join(PRD, 'surface', 'region_mapping_start_'+label+'.txt'), region_mapping_new, fmt='%d', newline=" ")


def correct_region_mapping(SUBJ_ID, PRD, label, REGION_MAPPING_CORR, DISPLAY, CHECK):

    print('automatically correct the new region mapping')

    texture = np.loadtxt(os.path.join(PRD, 'surface', 'region_mapping_start_'+label+'.txt'))[:20000]
    vert = np.loadtxt(os.path.join(PRD, SUBJ_ID, 'surface', 'vertices.txt'))
    trian = np.loadtxt(os.path.join(PRD, SUBJ_ID, 'surface', 'triangles.txt'))
    for ind in range(10):
        # print(ind)
        new_texture = deepcopy(texture)
        labels = np.unique(texture)
        for ilab in labels:
            iverts = (np.nonzero(texture==ilab)[0]).tolist()
            if len(iverts)>0:
                for inode in iverts:
                    iall = trian[np.nonzero(trian==inode)[0]].flatten().tolist()
                    ineig = np.unique(list(filter(lambda x: x!=inode, iall))).astype('int')
                    ivals = np.array(Counter(texture[ineig]).most_common()).astype('int')                
                    if ivals[np.nonzero(ivals[:,0]==ilab), 1].shape[1]==0:
                        new_texture[inode] = Counter(texture[ineig]).most_common(1)[0][0]
                    elif ivals[np.nonzero(ivals[:,0] == ilab), 1][0,0] < float(REGION_MAPPING_CORR) * len(ineig):
                        new_texture[inode] = Counter(texture[ineig]).most_common(1)[0][0]
        texture = deepcopy(new_texture) 

    np.savetxt(os.path.join(PRD, 'surface', 'region_mapping_not_corrected_'+label+'.txt'), new_texture)


def breadth_first_search(texture, iposi, itrian, ilab):
    queue = [iposi]
    V = [iposi]
    while len(queue) > 0:
        iQ = queue.pop()
        iedges = list(itrian[np.argwhere(itrian==iQ)[:,0]].flatten())
        while len(iedges)>0:
            ineig = iedges.pop()
            if ineig not in V and texture[ineig]==ilab:
                V.append(ineig)
                queue.append(ineig)
    return V


def calculate_connected(texture, vert, trian):
    "find if the regions are connected components using Breadth-first seach"
    labels = np.unique(texture) 
    res = []
    for ilab in labels:
        ipos = np.nonzero(texture==ilab)
        ivert = vert[ipos]
        itrian=[]
        for itri in np.nonzero(texture==ilab)[0].tolist():
            itrian.extend(trian[np.nonzero(trian==itri)[0]])
        itrian = np.array(itrian).astype('int')
        V = breadth_first_search(texture, ipos[0][0], itrian, ilab)
        res.append((ilab, ivert.shape[0]-len(V)))
    return res


def find_both_components(texture, vert, trian, ilab):
    " find the subgraphs"
    V = []
    V_done = []
    ipos = np.nonzero(texture==ilab)
    ivert = vert[ipos]
    itrian=[]
    for itri in ipos[0].tolist():
        itrian.extend(trian[np.nonzero(trian==itri)[0]])
    itrian = np.array(itrian).astype('int')
    # regions
    for istart in range(len(ipos[0])):
        if ipos[0][istart] not in V_done:
            V.append(breadth_first_search(texture, ipos[0][istart], itrian, ilab))
            V_done.extend(breadth_first_search(texture, ipos[0][istart], itrian, ilab))
    return V


def correct_sub_region(texture, trian, Vw):
    "correct the region mapping for the chosen component"
    new_texture = np.copy(texture)
    icount = 0
    while len(Vw)>0:
        iVw = Vw.pop()
        itrian = trian[np.nonzero(trian==iVw)[0]].flatten().astype('int').tolist()
        ir = list(filter(lambda x : new_texture[x] != new_texture[iVw], itrian))
        if len(ir)>0:
            new_texture[iVw] = new_texture[Counter(ir).most_common(1)[0][0]] 
        else:
            if icount<1000: 
                Vw.insert(0, iVw)
                icount +=1
            else:
                # TODO: good error message
                print('error in correction')
                import pdb; pdb.set_trace()
    return new_texture


def check_region_mapping(texture, vert, trian, ilab):
    "drawing the region"
    ipos = np.nonzero(texture==ilab)
    itrian=[]
    for itri in ipos[0].tolist():
        itrian.extend(trian[np.nonzero(trian==itri)[0]])
    itrian = np.array(itrian).astype('int')
    bool_itrian = np.in1d(itrian, ipos[0]).reshape(itrian.shape[0], 3)
    itrian[np.nonzero(bool_itrian == False)] = 0
    citri = vstack([vstack([itrian[:,0], itrian[:,1]]).T, vstack([itrian[:,1],itrian[:,2]]).T, vstack([itrian[:,2],itrian[:,0]]).T])
    bcitri = (citri!=0).sum(1)
    valp = citri[bcitri==2]
    fig = figure(figsize=(15, 15))
    fig.suptitle('region ' + str(int(ilab)))
    ax = fig.add_subplot(111, projection='3d') 
    xlabel('x')
    ylabel('y')
    for iv in np.arange(valp.shape[0]):
        ax.plot(vert[valp[iv], 0], vert[valp[iv], 1], vert[valp[iv], 2])
    # old function
    # xitrians = vert[np.hstack((itrian, itrian[:,0][:, newaxis])), 0]
    # yitrians = vert[np.hstack((itrian, itrian[:,0][:, newaxis])), 1]
    # zitrians = vert[np.hstack((itrian, itrian[:,0][:, newaxis])), 2]
    # fig = figure(figsize=(15, 15))
    # fig.suptitle('region ' + int(ilab))
    # ax = fig.add_subplot(111, projection='3d') 
    # for iv in range(xitrians.shape[0]): 
    #     ax.plot(xitrians[iv], yitrians[iv], zitrians[iv], alpha=0.4)
    # ax.scatter(vert[ipos, 0], vert[ipos, 1], vert[ipos, 2], c='b', s=45)
    show()


def main_check_region_mapping(subj_id, path_data, label, REGION_MAPPING_CORR,
	                          DISPLAY, CHECK):
        
    print('check and correct the new region mapping')
    vert = np.loadtxt(os.path.join(PRD, SUBJ_ID, 'surface', 'vertices.txt'))
    trian = np.loadtxt(os.path.join(PRD, SUBJ_ID, 'surface', 'triangles.txt'))
    texture = np.loadtxt(os.path.join(PRD, 'surface', 'region_mapping_not_corrected_'+label+'.txt'))[:20000]
    new_texture = deepcopy(texture)

    for i in range(3): # to remove all multiple components
        res = np.array(calculate_connected(new_texture, vert, trian))
        wrong_labels = res[np.where(res[:,1]>0.),0][0]
        if len(wrong_labels)==0:
            print('everything is fine, continuing')
        else:
            print("WARNING, some region have several components")
            print("ok if you only have the region 0.0 component after the final iteration")
            for iwrong in wrong_labels:
                # TODO: handle more than two components
                V = find_both_components(new_texture, vert, trian, iwrong)
                if iwrong==0.0:
                    choice_user = np.argsort([len(i) for i in V])[-2:]
                else:
                    choice_user = np.argsort([len(i) for i in V])[-1:]
                print(str(iwrong) + ' : ')
                for i in range(len(V)):
                    print(str(len(V[i])) + ' nodes')
                to_do = [i for i in range(len(V)) if i not in choice_user]
                for j in to_do:
                    new_texture  =  correct_sub_region(new_texture, trian, V[j])

    np.savetxt(os.path.join(PRD, SUBJ_ID, 'region_mapping_' + label +'.txt'), new_texture, fmt='%d', newline=' ')



if __name__ == '__main__':

    # config variables
    PRD = os.environ['PRD']
    CHECK = os.environ['CHECK']
    SUBJ_ID = os.environ['SUBJ_ID']
    REGION_MAPPING_CORR = os.environ['REGION_MAPPING_CORR']
    N_SUBREGIONS = os.environ['N_SUBREGIONS']
    N_SUBREGIONS = int(N_SUBREGIONS)
    if "DISPLAY" in os.environ:
        DISPLAY = os.environ['DISPLAY']
    else:
        DISPLAY = ""

    img = nib.load(os.path.join(PRD, 'connectivity', 'aparcaseg_2_diff_' + str(K) +'.nii.gz'))
    img_aff = nib.load(os.path.join(PRD, 'connectivity', 'aparc+aseg.nii.gz'))
    
    subcortical = [isub for isub in range(int(K)*70+1, int(K)*70+18)]

    surface_projection(SUBJ_ID, PRD, img, img_aff, N_SUBREGIONS, subcortical, DISPLAY, CHECK)
    
    correct_region_mapping(SUBJ_ID, PRD, N_SUBREGIONS, REGION_MAPPING_CORR, DISPLAY, CHECK)
    
    main_check_region_mapping(SUBJ_ID, PRD, N_SUBREGIONS, REGION_MAPPING_CORR, DISPLAY, CHECK)
        
