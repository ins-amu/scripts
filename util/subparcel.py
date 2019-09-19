import nibabel as nib
import os
import numpy as np
import copy
from sklearn import cluster
import scipy

# This script take a volumetric parcellation and divide each regions into subregions.
# The optimal number of subdivisions for each region is obtained so that we get the 
# smallest spread of the distribution of subregion volumes across all subregions after the subdivision

PRD = os.environ['PRD']
SUBJ_ID = os.environ['SUBJ_ID']
N_SUBREGIONS = os.environ['N_SUBREGIONS']
N_SUBREGIONS = int(N_SUBREGIONS)
PARCEL = os.environ['PARCEL']
CHECK = os.environ['CHECK']
if "DISPLAY" in os.environ:
    DISPLAY = os.environ['DISPLAY']
else:
    DISPLAY = ""

# check that N_SUBREGIONS is even
if N_SUBREGIONS%2==1:
    raise ValueError('N_SUBREGIONS must be an even number')


img = nib.load(os.path.join(PRD, 'connectivity', 'aparcaseg_2_diff.nii.gz'))
data_parcellation = img.get_data()

data_parcellation_cortical_only = copy.deepcopy(data_parcellation)
list_subcortical_regions = [2, 4, 5, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 
                            18, 24, 26, 28, 30, 31, 41, 43, 44, 46, 47, 49,
                            50, 51, 52, 53, 54, 58, 60, 62, 63, 72, 77, 80, 
                            85, 251, 252, 253, 254, 255]
for subcortical_region in list_subcortical_regions:
    data_parcellation_cortical_only[data_parcellation==subcortical_region] = 0

nib.save(img, os.path.join(PRD, 'connectivity', 'aparcaseg_2_diff_cortical_only.nii.gz'))

# check that N_SUBREGIONS is higher than the number of regions in the atlas
if N_SUBREGIONS<=np.unique(data_parcellation_cortical_only).shape[0]:
    raise ValueError('N_SUBREGIONS must higher than the number of regions in the original atlas')

# TODO subparcellation of subcortical structures
# TODO check that the whole cortex is in the img, otherwise add a border line of zeros; 
# that was necessary for the matlab version, probably not necessary anymore

## Find subparcellation that minimizes standard deviation of the distribution of region's sizes

# initialization
regions = np.unique(data_parcellation_cortical_only)[1:] # to remove label 0
# we want to subparcellate both hemisphere by the same number of regions
# TODO: improve for generic parcellation
if PARCEL=='desikan':
    regions_right = regions[(regions>=2000)]
elif PARCEL=='destrieux':
    regions_right = regions[(regions>=12000)]
else:
    raise ValueError('parcellation scheme not implemented, have a look in the code to add it')
l_idx_subdivisions = np.zeros(regions_right.shape[0], dtype=int)
l_nb_subdivisions = np.array([1, 2, 3, 4, 6, 8, 12, 16, 18, 24])
size_vol_parcellations = [data_parcellation_cortical_only[data_parcellation_cortical_only==region].shape[0] for region in regions_right]

# iteration right
# while(sum(l_nb_subdivisions)<N_SUBREGIONS/2):
#     curr_std = []
#     for iregion in np.arange(regions_right.shape[0]):
#         curr_size_vol_parcellation = []
#         for jregion in np.arange(regions_right.shape[0]):
#             if iregion!=jregion:
#                 curr_size_vol_parcellation.extend([size_vol_parcellations[jregion]/(l_nb_subdivisions[jregion]) for _ in np.arange(l_nb_subdivisions[jregion])])
#             elif iregion==jregion:
#                 curr_size_vol_parcellation.extend([size_vol_parcellations[iregion]/(l_nb_subdivisions[iregion]+1) for _ in np.arange(l_nb_subdivisions[iregion]+1)])
#         curr_std.append(np.std(curr_size_vol_parcellation))
#     l_nb_subdivisions[np.argmin(curr_std)] += 1
#     import pdb; pdb.set_trace()
# l_nb_subdivisions = np.tile(l_nb_subdivisions, (2,)) # for other hemisphere

while np.sum(l_nb_subdivisions[l_idx_subdivisions])<N_SUBREGIONS/2:
    curr_std = []
    for iregion in np.arange(regions_right.shape[0]):
        curr_size_vol_parcellation = []
        for jregion in np.arange(regions_right.shape[0]):
            if iregion!=jregion:
                curr_size_vol_parcellation.extend([size_vol_parcellations[jregion]/(l_nb_subdivisions[l_idx_subdivisions[jregion]]) for _ in np.arange(l_nb_subdivisions[l_idx_subdivisions[jregion]])])
            elif iregion==jregion:
                if l_idx_subdivisions[iregion]<9:
                    curr_size_vol_parcellation.extend([size_vol_parcellations[iregion]/(l_nb_subdivisions[l_idx_subdivisions[iregion]+1]) for _ in np.arange(l_nb_subdivisions[l_idx_subdivisions[iregion]+1])])
                else:
                    print('too many subregions, keep the same')
                    curr_size_vol_parcellation.extend([size_vol_parcellations[iregion]/(l_nb_subdivisions[l_idx_subdivisions[iregion]]) for _ in np.arange(l_nb_subdivisions[l_idx_subdivisions[iregion]])])
        curr_std.append(np.std(curr_size_vol_parcellation))
    iregion_update = 0
    while l_nb_subdivisions[l_idx_subdivisions[np.argsort(curr_std)[iregion_update]]+1]-l_nb_subdivisions[l_idx_subdivisions[np.argsort(curr_std)[iregion_update]]] > N_SUBREGIONS/2 - np.sum(l_nb_subdivisions[l_idx_subdivisions]):
        iregion_update +=1
    l_idx_subdivisions[np.argsort(curr_std)[iregion_update]] += 1
l_idx_subdivisions = np.tile(l_idx_subdivisions, (2,)) # for other hemisphere

if np.sum(l_nb_subdivisions[l_idx_subdivisions])!=N_SUBREGIONS:
    raise AssertionError('The number of subregions found by the algorihtm is different from the number of subregions required by the user')

## Finding same size subparcels for each region
new_data_parcellation = copy.deepcopy(data_parcellation_cortical_only)
for iregion in range(regions.shape[0]):
    l_subdivide_iterations = [[1], [2], [3], [2, 2], [2, 3], [2, 2, 2], [2, 2, 3], [2, 2, 2, 2], [2, 3, 3], [2, 2, 2, 3]]
    print('equalizing subparcels for region:' + str(iregion))

    l_curr_subregion = [regions[iregion]]
    for i_div in l_subdivide_iterations[l_idx_subdivisions[iregion]]:
        print(i_div)
        
        tmp_new_data_parcellation = copy.deepcopy(new_data_parcellation)
        for ii_curr_subregion, i_curr_subregion in enumerate(l_curr_subregion):
            # initialization: find k-mean centers, then iteratively grow the regions from it
            iregion_voxels = np.array(np.nonzero(new_data_parcellation==i_curr_subregion)).T
            k_means = cluster.MiniBatchKMeans(n_clusters=i_div)
            k_means.fit(np.array(iregion_voxels))
            cluster_centers = k_means.cluster_centers_
            n_clusters = cluster_centers.shape[0]
            cluster_labels = k_means.labels_
            _, subregion_counts = np.unique(k_means.labels_, return_counts=True)

            # TODO, be sure that the regions are in the same order across patients
            # iteration, reassign regions n the k-means to obtain same number of voxels in each subregion
            error_margin = np.sum(k_means.counts_)/(100*i_div) # 1%
            print(np.int(error_margin/2))
            while np.any((subregion_counts<np.mean(subregion_counts)-error_margin) | (subregion_counts>np.mean(subregion_counts)+error_margin)):
            #for i in range(100):
                print(subregion_counts)
                # estimate new cluster center
                for i_cluster in range(n_clusters):
                    cluster_centers[i_cluster] = np.mean(iregion_voxels[cluster_labels==i_cluster], 0) 
                # for biggest cluster, assign to closest region according to cost function
                # TODO: add a constraint on the distance to the center of the cluster
                i_smallest_clusters = np.where(subregion_counts<=(subregion_counts.min()+1))[0]
                i_bigger_clusters = np.nonzero(subregion_counts>(subregion_counts.min()+1))[0]
                candidate_clusters, candidate_points, candidate_new_cluster = [], [], []
                for i_smallest_cluster in i_smallest_clusters:
                    indices_smallest = np.where(cluster_labels==i_smallest_cluster)[0]
                    for i_bigger_cluster in i_bigger_clusters:
                        indices_bigger = np.where(cluster_labels==i_bigger_cluster)[0]
                        dist = scipy.spatial.distance.cdist(iregion_voxels[indices_smallest], iregion_voxels[indices_bigger])
                        i_cluster_min = np.amin(dist)
                        if i_cluster_min < 2:
                            candidate_clusters.append(i_bigger_cluster)
                            candidate_points.append(indices_bigger[np.where(dist==i_cluster_min)[1]])
                            candidate_new_cluster.append(i_smallest_cluster)
                changing_voxels_all = candidate_points[np.argmax(subregion_counts[candidate_clusters])]
                new_cluster = candidate_new_cluster[np.argmax(subregion_counts[candidate_clusters])]
                changing_voxels_closest = np.argsort(np.sum((iregion_voxels[changing_voxels_all]-cluster_centers[new_cluster])**2, 1))[:np.int(error_margin)] 
                #changing_region = np.random.choice(candidate_points[np.argmax(subregion_counts[candidate_clusters])], np.int(error_margin))
                cluster_labels[changing_voxels_all[changing_voxels_closest]] = new_cluster
                                
                # new subregion count
                _, subregion_counts = np.unique(cluster_labels, return_counts=True)
            
            tmp_new_data_parcellation[new_data_parcellation==i_curr_subregion] =  regions[iregion]*100 + ii_curr_subregion*10 + cluster_labels

            if CHECK=="yes":
                # figure test
                import matplotlib.pyplot as plt
                from mpl_toolkits.mplot3d import Axes3D
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
                ax.scatter(iregion_voxels[:,0], iregion_voxels[:,1], iregion_voxels[:,2], c=np.array(['r', 'b', 'g', 'y', 'k', 'm', 'c']*10)[cluster_labels.astype('int')], s=50, alpha=0.4)
                plt.show()

        new_data_parcellation =  copy.deepcopy(tmp_new_data_parcellation)
        l_curr_subregion = np.unique(new_data_parcellation[data_parcellation_cortical_only==regions[iregion]])

# check if loop didn't finish completely
if (np.unique(new_data_parcellation).shape[0]-1)!=N_SUBREGIONS:
    raise AssertionError('the number of subregions obtained by the algorithm is different from the number of subregions asked by the user')

# relabel regions by contiguous integer starting from 1
map_mat = np.vstack([np.unique(new_data_parcellation), np.arange(np.unique(new_data_parcellation).shape[0])])
data_subparcellation  = copy.deepcopy(new_data_parcellation)
for isub_reg, sub_reg in enumerate(np.unique(new_data_parcellation)):
    vox = np.nonzero(new_data_parcellation==sub_reg)
    data_subparcellation[vox] = isub_reg

last_reg = np.max(data_subparcellation)
# we only readd some subcortical regions
list_new_subcortical_regions = [16, 8, 10, 11, 12, 13, 17, 18, 26, 47, 49, 50, 51, 52, 53, 54, 58]
for isub_reg, subcortical_region in enumerate(list_new_subcortical_regions):
    data_subparcellation[data_parcellation==subcortical_region] = isub_reg + 1 + last_reg
new_hdr = img.header.copy()
new_hdr.set_data_dtype(np.uint8)
new_img = nib.Nifti1Image(data_subparcellation.astype(np.uint8), img.affine, new_hdr)
nib.save(new_img, os.path.join(PRD, 'connectivity', 'aparcaseg_2_diff_' + str(N_SUBREGIONS) + '.nii.gz'))    

## compute centers and orientations
ref_table = np.loadtxt('share/reference_table_' + PARCEL + '.csv', dtype='bytes', skiprows=1, delimiter=',').astype('str')

name_region = ref_table[:,1];
cortical = ref_table[:, 7].astype(bool);

name_region_subcortical = name_region[~cortical][1:] # we remove region unknown
name_region_cortical = name_region[cortical]
name_subregion_cortical, idx_subregion_cortical = [], []
for iregion, region in enumerate(name_region_cortical):
    name_subregion_cortical.extend([region]*l_nb_subdivisions[l_idx_subdivisions[iregion]])
    idx_subregion_cortical.extend([str(i) for i in range(l_nb_subdivisions[l_idx_subdivisions[iregion]])])

average_orientations = np.loadtxt(os.path.join(PRD, SUBJ_ID, 'connectivity', 'average_orientations.txt'));
average_orientations_subcortical = average_orientations[~cortical][1:]
average_orientations_cortical = average_orientations[cortical]
average_subregion_cortical = []
for iregion, region in enumerate(average_orientations_cortical):
    average_subregion_cortical.extend([region]*l_nb_subdivisions[l_idx_subdivisions[iregion]])

list_region = np.unique(data_subparcellation)[1:]
centres = np.zeros((list_region.shape[0], 4), dtype=object)
orientation_divided = np.zeros((list_region.shape[0], 3));

for iregion, region in enumerate(list_region[:list_region.shape[0] - len(list_new_subcortical_regions)]):
    pxl_iregion = np.nonzero(data_subparcellation==region)
    centres[iregion, 1:] = np.mean(pxl_iregion, 1)
    centres[iregion, 0] = name_subregion_cortical[iregion] + '_' + idx_subregion_cortical[iregion]
    orientation_divided[iregion] = average_subregion_cortical[iregion]

for iregion, region in enumerate(list_region[list_region.shape[0] - len(list_new_subcortical_regions):]): 
    pxl_iregion = np.nonzero(data_subparcellation==iregion); 
    centres[iregion + N_SUBREGIONS, 1:] = np.mean(pxl_iregion, 1);
    centres[iregion + N_SUBREGIONS, 0] = name_region_subcortical[iregion]; 
    orientation_divided[iregion + N_SUBREGIONS] = average_orientations_subcortical[iregion]; 


np.savetxt(os.path.join(PRD, SUBJ_ID, 'connectivity_' + str(N_SUBREGIONS), 'centres.txt'), centres, delimiter=" ", fmt="%s %.4f %.4f %.4f" )
np.savetxt(os.path.join(PRD, SUBJ_ID, 'connectivity_' + str(N_SUBREGIONS), 'average_orientations.txt'), orientation_divided)

# save corr_mat
np.savetxt(os.path.join(PRD, 'connectivity', 'corr_mat_' + str(N_SUBREGIONS) + '.txt'), np.vstack([list_region, list_region]).T)



    # # find the n
    # # TODO, be sure that the regions are in the same order across patients
    # # iteration, reassign regions n the k-means to obtain same number of voxels in each subregion
    # error_margin = np.sum(k_means.counts_)/(100*l_nb_subdivisions[iregion]) # 1%
    # while np.any((subregion_counts<np.mean(subregion_counts)-error_margin) | (subregion_counts>np.mean(subregion_counts)+error_margin)):
    # #for i in range(100):
    #     print(subregion_counts)
    #     # estimate new cluster center
    #     for i_cluster in range(n_clusters):
    #         cluster_centers[i_cluster] = np.mean(iregion_voxels[cluster_labels==i_cluster], 0) 

    #     # for biggest cluster, assign to closest region according to cost function
    #     i_smallest_clusters = np.where(subregion_counts<=(subregion_counts.min()+1))[0]
    #     i_bigger_clusters = np.nonzero(subregion_counts>(subregion_counts.min()+1))[0]
    #     candidate_clusters, candidate_points, candidate_new_cluster = [], [], []
    #     for i_smallest_cluster in i_smallest_clusters:
    #         indices_smallest = np.where(cluster_labels==i_smallest_cluster)[0]
    #         for i_bigger_cluster in i_bigger_clusters:
    #             indices_bigger = np.where(cluster_labels==i_bigger_cluster)[0]
    #             dist = scipy.spatial.distance.cdist(iregion_voxels[indices_smallest], iregion_voxels[indices_bigger])
    #             i_cluster_min = np.amin(dist)
    #             if i_cluster_min < 2:
    #                 candidate_clusters.append(i_bigger_cluster)
    #                 candidate_points.append(indices_bigger[np.where(dist==i_cluster_min)[1]])
    #                 candidate_new_cluster.append(i_smallest_cluster)
    #     changing_region = np.random.choice(candidate_points[np.argmax(subregion_counts[candidate_clusters])])
    #     cluster_labels[changing_region] = candidate_new_cluster[np.argmax(subregion_counts[candidate_clusters])]
        
    #     # new subregion count
    #     _, subregion_counts = np.unique(cluster_labels, return_counts=True)

#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     cm = plt.get_cmap('jet')
#     cs = range(iregion_voxels.shape[0])
#     cNorm = matplotlib.colors.Normalize(vmin=min(cs), vmax=max(cs))
#     scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
#     ax.scatter(iregion_voxels[:,0], iregion_voxels[:,1], iregion_voxels[:,2], s=50, alpha=0.1)
#     ax.scatter(iregion_voxels[center,0], iregion_voxels[center,1], iregion_voxels[center,2], c='r', s=200, alpha=1)
#     #ax.scatter(iregion_voxels[:,0], iregion_voxels[:,1], iregion_voxels[:,2], c=scalarMap.to_rgba(cs), s=50, alpha=0.1)

#     figure()
#     scatter(Y[:,0], Y[:,1], c=scalarMap.to_rgba(cs))
#     scalarMap.set_array(cs)
#     show()

#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     ax.scatter(Y[:,0], Y[:,1], Y[:,2], s=50, alpha=0.1)
#     plt.show()

#     # figure test
#     import matplotlib.pyplot as plt
#     from mpl_toolkits.mplot3d import Axes3D
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     ax.scatter(iregion_voxels[:,0], iregion_voxels[:,1], iregion_voxels[:,2], c=np.array(['r', 'b', 'g', 'y', 'k']*10)[cluster_labels], s=50, alpha=0.1)
#     #ax.scatter(iregion_voxels[indices_bigger,0], iregion_voxels[indices_bigger,1], iregion_voxels[indices_bigger,2], c=np.array(['r', 'b', 'g', 'y', 'k']*10)[cluster_labels[indices_bigger]], s=10, alpha=0.1)
#     #ax.scatter(iregion_voxels[indices_bigger[changing_region],0], iregion_voxels[indices_bigger[changing_region],1], iregion_voxels[indices_bigger[changing_region],2], c=np.array(['r']), s=100)
#     ax.scatter(cluster_centers[:,0], cluster_centers[:,1], cluster_centers[:,2], c= np.array(['k']), s=100)
#     plt.show()


# for iregion in range(10):#regions.shape[0]):
#     print('equalizing subparcels for region:' + str(iregion))
#     # initialization: find most separated points from locally linear embedding

#     iregion_voxels = np.array(np.nonzero(data_parcellation==regions[iregion])).T
#     k_means = cluster.MiniBatchKMeans(n_clusters=l_nb_subdivisions[iregion])
#     k_means.fit(np.array(iregion_voxels))
#     cluster_centers = k_means.cluster_centers_
#     n_clusters = cluster_centers.shape[0]
#     cluster_labels = k_means.labels_
#     _, subregion_counts = np.unique(k_means.labels_, return_counts=True)


#     # TODO, be sure that the regions are in the same order across patients
#     # iteration, reassign regions n the k-means to obtain same number of voxels in each subregion
#     error_margin = np.sum(k_means.counts_)/(100*l_nb_subdivisions[iregion]) # 1%
#     while np.any((subregion_counts<np.mean(subregion_counts)-error_margin) | (subregion_counts>np.mean(subregion_counts)+error_margin)):
#     #for i in range(100):
#         print(subregion_counts)
#         # estimate new cluster center
#         for i_cluster in range(n_clusters):
#             cluster_centers[i_cluster] = np.mean(iregion_voxels[cluster_labels==i_cluster], 0) 

#         # for biggest cluster, assign to closest region according to cost function
#         i_smallest_clusters = np.where(subregion_counts<=(subregion_counts.min()+1))[0]
#         i_bigger_clusters = np.nonzero(subregion_counts>(subregion_counts.min()+1))[0]
#         candidate_clusters, candidate_points, candidate_new_cluster = [], [], []
#         for i_smallest_cluster in i_smallest_clusters:
#             indices_smallest = np.where(cluster_labels==i_smallest_cluster)[0]
#             for i_bigger_cluster in i_bigger_clusters:
#                 indices_bigger = np.where(cluster_labels==i_bigger_cluster)[0]
#                 dist = scipy.spatial.distance.cdist(iregion_voxels[indices_smallest], iregion_voxels[indices_bigger])
#                 i_cluster_min = np.amin(dist)
#                 if i_cluster_min < 2:
#                     candidate_clusters.append(i_bigger_cluster)
#                     candidate_points.append(indices_bigger[np.where(dist==i_cluster_min)[1]])
#                     candidate_new_cluster.append(i_smallest_cluster)
#         changing_region = np.random.choice(candidate_points[np.argmax(subregion_counts[candidate_clusters])])
#         cluster_labels[changing_region] = candidate_new_cluster[np.argmax(subregion_counts[candidate_clusters])]
        
#         # new subregion count
#         _, subregion_counts = np.unique(cluster_labels, return_counts=True)


#     # figure test
#     import matplotlib.pyplot as plt
#     from mpl_toolkits.mplot3d import Axes3D
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     ax.scatter(iregion_voxels[:,0], iregion_voxels[:,1], iregion_voxels[:,2], c=np.array(['r', 'b', 'g', 'y', 'k', 'm', 'c']*10)[cluster_labels.astype('int')], s=50, alpha=0.4)
#     #ax.scatter(iregion_voxels[indices_bigger,0], iregion_voxels[indices_bigger,1], iregion_voxels[indices_bigger,2], c=np.array(['r', 'b', 'g', 'y', 'k']*10)[cluster_labels[indices_bigger]], s=10, alpha=0.1)
#     #ax.scatter(iregion_voxels[indices_bigger[changing_region],0], iregion_voxels[indices_bigger[changing_region],1], iregion_voxels[indices_bigger[changing_region],2], c=np.array(['r']), s=100)
#     #ax.scatter(cluster_centers[:,0], cluster_centers[:,1], cluster_centers[:,2], c= np.array(['k']), s=100)
#     plt.show()



# # TODO surface parcellation

# ## Finding same size subparcels for each region
# for iregion in range(10):#regions.shape[0]):
#     l_subdivide_iterations = [[1], [2], [3], [2, 2], [2, 3], [2, 2, 2], [2, 2, 3], [2, 2, 2, 2], [2, 3, 3], [2, 2, 2, 3]]
#     print('equalizing subparcels for region:' + str(iregion))

#     new_data_parcellation= copy.deepcopy(data_parcellation)
#     l_curr_subregion = [regions[iregion]]
#     for i_div in l_subdivide_iterations[l_idx_subdivisions[iregion]]:

#         for i_curr_subregion in l_curr_subregion:
#             # initialization: find k-mean centers, then iteratively grow the regions from it
#             iregion_voxels = np.array(np.nonzero(new_data_parcellation==i_curr_subregion)).T
#             clusters = cluster.SpectralClustering(n_clusters=i_div)
#             clusters.fit(np.array(iregion_voxels))
#             cluster_centers = [np.mean(iregion_voxels[clusters.labels_==icluster], axis=0) for icluster in range(l_nb_subdivisions[iregion])]
#             n_clusters = np.array(cluster_centers).shape[0]


#             target_size = np.floor(iregion_voxels.shape[0]/i_div)
#             ltarget_size = np.ones(l_nb_subdivisions[iregion]) * target_size
#             remaining_voxels = iregion_voxels.shape[0]%l_nb_subdivisions[iregion]
#             ltarget_size[0:remaining_voxels] += 1

#             cluster_labels = np.zeros(iregion_voxels.shape[0])
#             curr_iregion_voxels = copy.deepcopy(iregion_voxels)
#             l_curr_subregion = []
#             for icluster in range(i_div):
#                 # compute distance of all remaining points to the center
#                 ldist = scipy.spatial.distance.pdist(curr_iregion_voxels)
#                 dist = scipy.spatial.distance.squareform(ldist)
#                 center = np.argmax(np.sum(dist, 1))
#                 cluster_labels[center] = icluster
#                 cluster_indices = np.argsort(dist[center])[1:ltarget_size[icluster]]
#                 cluster_labels[cluster_indices] = icluster
#                 print(cluster_indices.shape)
#                 #import pdb; pdb.set_trace()
#                 curr_iregion_voxels = np.delete(curr_iregion_voxels, np.hstack([center, cluster_indices]), axis=0)
#                 print(curr_iregion_voxels.shape)
#                 l_curr_subregion.append(region[iregion]*100 + icluster)
            
#             new_data_parcellation[new_data_parcellation==regions[iregion]] = cluster_labels

#             # figure test
#             import matplotlib.pyplot as plt
#             from mpl_toolkits.mplot3d import Axes3D
#             fig = plt.figure()
#             ax = fig.add_subplot(111, projection='3d')
#             ax.scatter(iregion_voxels[:,0], iregion_voxels[:,1], iregion_voxels[:,2], c=np.array(['r', 'b', 'g', 'y', 'k', 'm', 'c']*10)[cluster_labels.astype('int')], s=50, alpha=0.4)
#             plt.show()
