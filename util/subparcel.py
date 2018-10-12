import nibabel as nib
import os
import numpy as np
import copy
from sklearn import cluster
import scipy

PRD = os.environ['PRD']
SUBJ_ID = os.environ['SUBJ_ID']
n_subregions = os.environ['N_SUBREGIONS']
CHECK = os.environ['CHECK']
if "DISPLAY" in os.environ:
    DISPLAY = os.environ['DISPLAY']
else:
    DISPLAY = ""


data_parcellation = nib.load(os.path.join(PRD, 'connectivity', 'aparcaseg_2_diff.nii.gz')).get_data()

list_subcortical_regions = [2, 4, 5, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 
                            18, 24, 26, 28, 30, 31, 41, 43, 44, 46, 47, 49,
                            50, 51, 52, 53, 54, 58, 60, 62, 63, 72, 77, 80, 
                            85, 251, 252, 253, 254, 255]
for subcortical_region in list_subcortical_regions:
    data_parcellation[data_parcellation==subcortical_region] = 0

# TODO subparcellation of subcortical structures
# TODO test if n_regions is list of numbers, if so skip next step
# TODO require that n_subregions is even

## Find subparcellation that minimizes standard deviation

# initialization
regions = np.unique(data_parcellation)[1:] # to remove label 0
regions_right = regions[(regions>=2000)] # we want to subparcellate both hemisphere by the same number of regions
l_nb_subdivisions = np.ones(regions_right.shape[0], dtype=int)
size_vol_parcellations = [data_parcellation[data_parcellation==region].shape[0] for region in regions_right]

# iteration right
while(sum(l_nb_subdivisions)<n_subregions/2):
    curr_std = []
    for iregion in np.arange(regions_right.shape[0]):
        curr_size_vol_parcellation = []
        for jregion in np.arange(regions_right.shape[0]):
            if iregion!=jregion:
                curr_size_vol_parcellation.extend([size_vol_parcellations[jregion]/(l_nb_subdivisions[jregion]) for _ in np.arange(l_nb_subdivisions[jregion])])
            elif iregion==jregion:
                curr_size_vol_parcellation.extend([size_vol_parcellations[iregion]/(l_nb_subdivisions[iregion]+1) for _ in np.arange(l_nb_subdivisions[iregion]+1)])
        curr_std.append(np.std(curr_size_vol_parcellation))
    l_nb_subdivisions[np.argmin(curr_std)] += 1
l_nb_subdivisions = np.tile(l_nb_subdivisions, (2,)) # for other hemisphere


## Finding same size subparcels for each region

for iregion in range(10):#regions.shape[0]):
    print('equalizing subparcels for region:' + str(iregion))
    # initialization: k-mean
    iregion_voxels = np.array(np.nonzero(data_parcellation==regions[iregion])).T
    k_means = cluster.MiniBatchKMeans(n_clusters=l_nb_subdivisions[iregion])
    k_means.fit(np.array(iregion_voxels))
    cluster_centers = k_means.cluster_centers_
    n_clusters = cluster_centers.shape[0]
    cluster_labels = k_means.labels_
    _, subregion_counts = np.unique(k_means.labels_, return_counts=True)


    # TODO, be sure that the regions are in the same order across patients
    # iteration, reassign regions n the k-means to obtain same number of voxels in each subregion
    error_margin = np.sum(k_means.counts_)/(100*l_nb_subdivisions[iregion]) # 1%
    while np.any((subregion_counts<np.mean(subregion_counts)-error_margin) | (subregion_counts>np.mean(subregion_counts)+error_margin)):
    #for i in range(100):
        print(subregion_counts)
        # estimate new cluster center
        for i_cluster in range(n_clusters):
            cluster_centers[i_cluster] = np.mean(iregion_voxels[cluster_labels==i_cluster], 0) 

        # for biggest cluster, assign to closest region according to cost function
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
        changing_region = np.random.choice(candidate_points[np.argmax(subregion_counts[candidate_clusters])])
        cluster_labels[changing_region] = candidate_new_cluster[np.argmax(subregion_counts[candidate_clusters])]
        
        # new subregion count
        _, subregion_counts = np.unique(cluster_labels, return_counts=True)


    # figure test
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(iregion_voxels[:,0], iregion_voxels[:,1], iregion_voxels[:,2], c=np.array(['r', 'b', 'g', 'y', 'k']*10)[cluster_labels], s=50, alpha=0.1)
    #ax.scatter(iregion_voxels[indices_bigger,0], iregion_voxels[indices_bigger,1], iregion_voxels[indices_bigger,2], c=np.array(['r', 'b', 'g', 'y', 'k']*10)[cluster_labels[indices_bigger]], s=10, alpha=0.1)
    #ax.scatter(iregion_voxels[indices_bigger[changing_region],0], iregion_voxels[indices_bigger[changing_region],1], iregion_voxels[indices_bigger[changing_region],2], c=np.array(['r']), s=100)
    ax.scatter(cluster_centers[:,0], cluster_centers[:,1], cluster_centers[:,2], c= np.array(['k']), s=100)
    plt.show()



# TODO surface parcellation