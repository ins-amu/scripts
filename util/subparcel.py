import nibabel as nib
import os
import numpy as np
import copy

PRD = os.environ['PRD']
SUBJ_ID = os.environ['SUBJ_ID']
n_subregions = os.environ['N_SUBREGIONS']

data_parcellation = nib.load(os.path.join(PRD, 'connectivity', 'aparcaseg_2_diff.nii.gz')).get_data()

list_subcortical_regions = [2, 4, 5, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 
                            18, 24, 26, 28, 30, 31, 41, 43, 44, 46, 47, 49,
                            50, 51, 52, 53, 54, 58, 60, 62, 63, 72, 77, 80, 
                            85, 251, 252, 253, 254, 255]
for subcortical_region in list_subcortical_regions:
    data_parcellations[data_parcellations==subcortical_region] = 0

# TODO subparcellation of subcortical structures
# TODO test if n_regions is list of numbers, if so skip next step
# Find subparcellation that minimizes standard deviation


# initialization
regions = np.unique(data_parcellations)[1:] # to remove label 0
regions_right = regions[(regions>=2000)] # TODO check - we want subparcellate by the same number of regions
l_nb_subdivisions = np.ones(regions_right.shape[0], dtype=int)
size_vol_parcellations = [data_parcellations[data_parcellations==region].shape[0] for region in regions_right]

# iteration right
while(sum(l_nb_subdivisions[1:region])<n_subregions/2):
    curr_std = []
    for iregion in np.arange(regions_right.shape[0]):
        curr_size_vol_parcellation = []
        for jregion in np.arange(regions.shape[0]):
            if iregion!=jregion:
                curr_size_vol_parcellation.extend([size_vol_parcellations[jregion]/(l_nb_subdivisions[jregion]) for _ in np.arange(l_nb_subdivisions[jregion])])
            elif iregion==jregion:
                curr_size_vol_parcellation.extend([size_vol_parcellations[iregion]/(l_nb_subdivisions[iregion]+1) for _ in np.arange(l_nb_subdivisions[iregion]+1)])
        curr_std.append(std(curr_size_vol_parcellation))
    l_nb_subdivisions[argmin(curr_std)] += 1
l_nb_subdivisions.extend(l_nb_subdivisions) # for otr hemisphere


# Finding same size subparcels for each region

for iregion in np.arange(regions.shape[0]):
    # initialization: k-mean
    iregion_voxels = np.nonzero(data_parcellation==regions[iregion])
    k_means = cluster.MiniBatchKMeans(n_clusters=l_nb_subdivisions[iregion])
    k_means.fit(np.array(iregion_voxels).T)
    cluster_centers = k_means.cluster_centers_
    n_clusters = cluster_centers.shape[0]
    cluster_labels = k_means.labels_

    # TODO, be sure that the regions are in the same order across patients
    # iteration, reassign regions n the k-means to obtain same number of voxels in each subregion
    cost_function = np.zeros((iregion_voxels[0].shape[0], n_clusters))

    error_margin = np.sum(k_means.count_)/(100*l_nb_subdivisions[iregion]) # 1%
    subregion_counts = k_means.counts_
    while np.any(subregion_counts<subregion_counts-error_margin | subregion_counts>subregion_counts+error_margin):
        # estimate new cluster center
        for i_cluster in range(n_cluster):
            cluster_centers[i_cluster] = np.mean(iregion_voxels[cluster_labels==icluster]) 
        # for biggest cluster, assign to closest region according to cost function
        i_biggest_cluster = argmax(subregion_counts)
        i_smaller_clusters = np.nonzero(subregions_count<subregions_count[i_biggest_cluster])
        distance_to_biggest_cluster = np.linalg.norm(np.array(iregion_voxels).T - cluster_centers[i_biggest_cluster], axis=1)
        # for the point that minimize the summed distance between the current cluster
        # and the closest second cluster
        for icluster in range(i_smaller_clusters):
             cost_function[i_smaller_clusters[icluster]] = np.linalg.norm(np.array(iregion_voxels[cluster_labels==icluster]).T - cluster_centers[icluster]) + distance_to_biggest_cluster
        [new_cluster, changing_region] = np.where(cost_function[i_smaller_clusters]==min(cost_function[i_smaller_clusters]), axis=1)
        cluster_labels[changing_region] = new_cluster



# TODO surface parcellation