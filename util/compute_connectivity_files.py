import numpy as np
import os

PRD = os.environ['PRD']
SUBJ_ID = os.environ['SUBJ_ID']
PARCEL = os.environ['PARCEL']

def compute_triangle_areas(vertices, triangles):
    """Calculates the area of triangles making up a surface."""
    tri_u = vertices[triangles[:, 1], :] - vertices[triangles[:, 0], :]
    tri_v = vertices[triangles[:, 2], :] - vertices[triangles[:, 0], :]
    tri_norm = np.cross(tri_u, tri_v)
    triangle_areas = np.sqrt(np.sum(tri_norm ** 2, axis=1)) / 2.0
    triangle_areas = triangle_areas[:, np.newaxis]
    return triangle_areas


def compute_region_areas(triangles_areas, vertex_triangles):
    avt = np.array(vertex_triangles)
    #NOTE: Slightly overestimates as it counts overlapping border triangles,
    #      but, not really a problem provided triangle-size << region-size.
    regs = map(set, avt)
    region_triangles = set.union(*regs)
    region_surface_area = triangle_areas[list(region_triangles)].sum()
    return region_surface_area


def compute_region_orientation(vertex_normals):
    average_orientation = np.zeros((1, 3))
    # Average orientation of the region
    orient = vertex_normals[:, :]
    avg_orient = np.mean(orient, axis=0)
    average_orientation = avg_orient / np.sqrt(np.sum(avg_orient ** 2))
    region_orientation = average_orientation
    return region_orientation


def compute_vertex_triangles(number_of_vertices, number_of_triangles, triangles):
    vertex_triangles = [[] for _ in range(number_of_vertices)]
    for k in range(number_of_triangles):
        vertex_triangles[triangles[k, 0]].append(k)
        vertex_triangles[triangles[k, 1]].append(k)
        vertex_triangles[triangles[k, 2]].append(k)
    return vertex_triangles


def compute_vertex_normals(number_of_vertices, vertex_triangles, triangles,
                           triangle_angles, triangle_normals, vertices):
    """
    Estimates vertex normals, based on triangle normals weighted by the
    angle they subtend at each vertex...
    """
    vert_norms = np.zeros((number_of_vertices, 3))
    bad_normal_count = 0
    for k in range(number_of_vertices):
        try:
            tri_list = list(vertex_triangles[k])
            angle_mask = triangles[tri_list, :] == k
            angles = triangle_angles[tri_list, :]
            angles = angles[angle_mask][:, np.newaxis]
            angle_scaling = angles / np.sum(angles, axis=0)
            vert_norms[k, :] = np.mean(angle_scaling * triangle_normals[tri_list, :], axis=0)
            # Scale by angle subtended.
            vert_norms[k, :] = vert_norms[k, :] / np.sqrt(np.sum(vert_norms[k, :] ** 2, axis=0))
            # Normalise to unit vectors.
        except (ValueError, FloatingPointError):
            # If normals are bad, default to position vector
            # A nicer solution would be to detect degenerate triangles and ignore their
            # contribution to the vertex normal
            vert_norms[k, :] = vertices[k] / np.sqrt(vertices[k].dot(vertices[k]))
            bad_normal_count += 1
    if bad_normal_count:
        print(" %d vertices have bad normals" % bad_normal_count)
    return vert_norms


def compute_triangle_angles(vertices, number_of_triangles, triangles):
    """
    Calculates the inner angles of all the triangles which make up a surface
    """
    verts = vertices
    # TODO: Should be possible with arrays, ie not nested loops...
    # A short profile indicates this function takes 95% of the time to compute normals
    # (this was a direct translation of some old matlab code)
    angles = np.zeros((number_of_triangles, 3))
    for tt in range(number_of_triangles):
        triangle = triangles[tt, :]
        for ta in range(3):
            ang = np.roll(triangle, -ta)
            angles[tt, ta] = np.arccos(np.dot(
                (verts[ang[1], :] - verts[ang[0], :]) /
                np.sqrt(np.sum((verts[ang[1], :] - verts[ang[0], :]) ** 2, axis=0)),
                (verts[ang[2], :] - verts[ang[0], :]) /
                np.sqrt(np.sum((verts[ang[2], :] - verts[ang[0], :]) ** 2, axis=0))))
    return angles


def compute_triangle_normals(triangles, vertices):
    """Calculates triangle normals."""
    tri_u = vertices[triangles[:, 1], :] - vertices[triangles[:, 0], :]
    tri_v = vertices[triangles[:, 2], :] - vertices[triangles[:, 0], :]
    tri_norm = np.cross(tri_u, tri_v)

    try:
        triangle_normals = tri_norm / np.sqrt(np.sum(tri_norm ** 2, axis=1))[:, np.newaxis]
    except FloatingPointError:
        #TODO: NaN generation would stop execution, however for normals this case could maybe be
        # handled in a better way.
        triangle_normals = tri_norm
    return triangle_normals


def compute_region_areas_cortex(triangle_areas, vertex_triangles, region_mapping, list_name):
    regions = np.unique(region_mapping)
    region_surface_area = np.zeros((list_name.shape[0], 1))
    avt = np.array(vertex_triangles)
    #NOTE: Slightly overestimates as it counts overlapping border triangles,
    #      but, not really a problem provided triangle-size << region-size.
    for k in regions:
        regs = map(set, avt[region_mapping == k])
        region_triangles = set.union(*regs)
        region_surface_area[k] = triangle_areas[list(region_triangles)].sum()
    return region_surface_area


def compute_region_orientation_cortex(vertex_normals, region_mapping, list_name):
    regions = np.unique(region_mapping)
    average_orientation = np.zeros((list_name.shape[0], 3))
    #Average orientation of the region
    for k in regions:
        orient = vertex_normals[region_mapping == k, :]
        avg_orient = np.mean(orient, axis=0)
        average_orientation[k, :] = avg_orient / np.sqrt(np.sum(avg_orient ** 2))
    return average_orientation


def compute_region_center_cortex(vertices, region_mapping, list_name):
    regions = np.unique(region_mapping)
    region_center= np.zeros((list_name.shape[0], 3))
    #Average orientation of the region
    for k in regions:
        vert = vertices[region_mapping == k, :]
        region_center[k, :] = np.mean(vert, axis=0)
    return region_center


if __name__ == '__main__':

    # Cortex
    # import data
    verts = np.loadtxt(os.path.join(PRD, SUBJ_ID, 'surface', 'vertices.txt'))
    tri = np.loadtxt(os.path.join(PRD, SUBJ_ID, 'surface', 'triangles.txt'))
    tri = tri.astype(int)
    region_mapping = np.loadtxt(os.path.join(PRD, SUBJ_ID, 'surface', 'region_mapping.txt')).astype(int)
    # save connectivity and tract length matrices
    weights = np.loadtxt(os.path.join(PRD, 'connectivity', 'weights.csv'))
    tract_lengths = np.loadtxt(os.path.join(PRD, 'connectivity', 'tract_lengths.csv'))
    weights = weights + weights.transpose() - np.diag(np.diag(weights))
    # add the first region
    weights = np.vstack([np.zeros((1, weights.shape[0])), weights])
    weights = np.hstack([np.zeros((weights.shape[0], 1)), weights])
    tract_lengths = tract_lengths + tract_lengths.transpose() # because diagonal nul 
    tract_lengths = np.vstack([np.zeros((1, tract_lengths.shape[0])), tract_lengths]) 
    tract_lengths = np.hstack([np.zeros((tract_lengths.shape[0], 1)), tract_lengths])
    np.savetxt(os.path.join(PRD, SUBJ_ID, 'connectivity', 'weights.txt'), weights, fmt='%d')
    np.savetxt(os.path.join(PRD, SUBJ_ID, 'connectivity', 'tract_lengths.txt'), tract_lengths, fmt='%.3f')

    # name of the centers
    list_name = np.loadtxt(open(os.path.join('share', 'reference_table_' + PARCEL + ".csv"), "r"), delimiter=",", skiprows=1, usecols=(1, ), dtype='str')

    # compute centers
    centers = compute_region_center_cortex(verts, region_mapping, list_name)

    # calculate average orientations
    number_of_vertices = int(verts.shape[0])
    number_of_triangles = int(tri.shape[0])
    vertex_triangles = compute_vertex_triangles(number_of_vertices, number_of_triangles,
                                        tri)
    triangle_normals = compute_triangle_normals(tri, verts)
    triangle_angles = compute_triangle_angles(verts, number_of_triangles, tri)
    vertex_normals = compute_vertex_normals(number_of_vertices, vertex_triangles,
                                            tri, triangle_angles,
                                            triangle_normals, verts)
    orientations = compute_region_orientation_cortex(vertex_normals, region_mapping, list_name)

    # compute areas
    triangle_areas = compute_triangle_areas(verts, tri)
    areas = compute_region_areas_cortex(triangle_areas, vertex_triangles, region_mapping, list_name)

    # subcorticals
    corr_table = np.loadtxt(open(os.path.join('share', 'reference_table_' + PARCEL + ".csv"), "r"), delimiter=",", skiprows=1, usecols=(0, 5))
    for val in ['16', '08', '10', '11', '12', '13', '17', '18', '26', '47', '49',
                '50', '51', '52', '53', '54', '58']:
        verts = np.loadtxt(os.path.join(PRD, 'surface', 'subcortical',
                                        'aseg_0'+str(val)+'_vert.txt'))
        tri = np.loadtxt(os.path.join(PRD, 'surface', 'subcortical',
                                        'aseg_0'+str(val)+'_tri.txt'))
        tri = tri.astype(int)

        curr_center = np.mean(verts, axis=0)
        indx = int(corr_table[np.nonzero(corr_table[:, 0] == np.int(val)), 1])
        centers[indx, :] = curr_center
        # Now calculate average orientations
        number_of_vertices = int(verts.shape[0])
        number_of_triangles = int(tri.shape[0])
        vertex_triangles = compute_vertex_triangles(number_of_vertices, number_of_triangles,
                                            tri)
        triangle_normals = compute_triangle_normals(tri, verts)
        triangle_angles = compute_triangle_angles(verts, number_of_triangles, tri)
        vertex_normals = compute_vertex_normals(number_of_vertices, vertex_triangles,
                                                tri, triangle_angles,
                                                triangle_normals, verts)
        average_orientation = compute_region_orientation(vertex_normals)
        orientations[indx, :] = average_orientation

        triangle_areas = compute_triangle_areas(verts, tri)
        region_areas = compute_region_areas(triangle_areas, vertex_triangles)
        areas[indx] = region_areas

    # save orientations and areas
    np.savetxt(os.path.join(PRD, SUBJ_ID, 'connectivity/areas.txt'), areas, fmt='%.2f')
    np.savetxt(os.path.join(PRD, SUBJ_ID, 'connectivity/average_orientations.txt'),
            orientations, fmt='%.2f %.2f %.2f')


    f = open(os.path.join(PRD, SUBJ_ID, 'connectivity/centres.txt'), 'w')
    for i, name in enumerate(list_name):
        f.write(str(name) +' ')
        for j in range(3):
            f.write('{:.4f} '.format(centers[i, j]))
        f.write('\n')
    f.close()

    # save cortical
    ref_table = np.loadtxt(open(os.path.join('share', 'reference_table_' + PARCEL + ".csv"), "rb"), delimiter=",", skiprows=1, usecols=(7,))
    np.savetxt(os.path.join(PRD, SUBJ_ID, 'connectivity/cortical.txt'), ref_table, fmt='%d')