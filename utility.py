
def extract_high(surface, rl):
    import numpy as np
    import os
    name_file = rl + surface
    with open(surface, 'r') as f:
        f.readline()
        nb_vert = f.readline().split(' ')[0]
        read_data = [[np.double(line.rstrip('\n').split()[0]),
                     np.double(line.rstrip('\n').split()[1]),
                     np.double(line.rstrip('\n').split()[2])] for line in f]

    a = np.array(read_data)
    vert_high = a[0:int(nb_vert), 0:3]
    tri_high = a[int(nb_vert):, 0:3]
    np.savetxt(rl + '_vertices_high.txt', vert_high, fmt='%.6f %.6f %.6f')
    tri_high = a[int(nb_vert):, 0:3]
    np.savetxt(rl +'_triangles_high.txt', tri_high, fmt='%d %d %d')
    return (map(os.path.abspath, [rl + '_vertices_high.txt', rl + '_triangles_high.txt'])

def txt2off(vertices, triangles, rl):
    import numpy as np
    import os

    vert = np.loadtxt(vertices)
    tri = np.loadtxt(triangles)

    with open(rl + '_high.off', 'w') as f: f.write('OFF\n') 
        f.write('{} {} {}\n'.format(int(vert.shape[0]), int(tri.shape[0]), 0)) 

    with open(rl + '_high.off', 'a') as f:
        np.savetxt(f, vert, fmt='%.6f')
        np.savetxt(f, np.hstack([np.ones((tri.shape[0],1))*3, tri]), fmt='%d')
     
    return os.path.abspath(rl + '_high.off')

def off2txt(surface, rl):
    import numpy as np
    import os

    surface
    with open(surface) as f:
        f.readline()
        num = f.readline().split(' ')

    vert = np.loadtxt(surface, skiprows=2, usecols=(0,1,2))  
    vert = vert[:int(num[0]), :]
    tri = np.loadtxt(surface, skiprows=int(num[0])+2, usecols=(1,2,3))  

    np.savetxt(rl + '_vertices_low.off', vert, fmt='%.4f')
    np.savetxt(rl + '_triangles_low.off', tri, fmt='%d')
    return (map(os.path.abspath, [rl + '_vertices_high.txt', rl + '_triangles_high.txt'])

def correct_region_mapping(region_mapping_not_corrected, vertices, triangles, rl,
                           region_mapping_corr=0.42):
    import os
    import sys
    region_mapping_corr = float(os.environ['region_mapping_corr'])
    from copy import deepcopy
    from pylab import *
    from collections import Counter

    texture = loadtxt(region_mapping_not_corrected)
    vert = loadtxt(vertices)
    trian = loadtxt(triangles)
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
    return os.path.abspath(rl + '_region_mapping_low.txt')

def reunify_both_regions(region_mapping_list, vertices_list, triangles_list, rl):
    import os
    from pylab import *
    lh_reg_map = np.loadtxt(region_mapping_list([0])
    lh_vert = np.loadtxt(vertices_list[0])
    lh_trian = np.loadtxt(triangles_list[0])
    rh_reg_map = np.loadtxt(region_mapping_list[1])
    rh_vert = np.loadtxt(vertices_list[1])
    rh_trian = np.loadtxt(triangles_list[1])
    vertices = vstack([lh_vert, rh_vert])
    triangles = vstack([lh_trian,  rh_trian + lh_vert.shape[0]])
    region_mapping = hstack([lh_reg_map, rh_reg_map])
    np.savetxt('region_mapping.txt', region_mapping, fmt='%d', newline=" ")
    np.savetxt('vertices.txt', vertices, fmt='%.2f')
    np.savetxt('triangles.txt', triangles, fmt='%d %d %d')
    return (map(os.path.abspath(), ['region_mapping.txt', 'vertices.txt', 'triangles.txt']))


class Aseg2SrfInputSpec(CommandLineSpec):
    in_subject_id = File(desc = "Subject FreeSurfer Id",
            argstr = '%d',
            exists = True,
            mandatory = True)


class Aseg2SrfOutputSpec(TraitedSpec):
    out_subcortical_surf_list = File(desc = "Output subcortical surfaces", exists = True)


class Aseg2Srf(CommandLine):
    input_spec = Aseg2SrfInputSpec
    output_spec = Aseg2SrfOutputSpec
    _cmd = './aseg2srf' 

    def _gen_subjects_dir(self):
        return os.getcwd()

    def _list_outputs(self):
        if isdefined(self.inputs.subjects_dir):
            subjects_dir = self.inputs.subjects_dir
        else:
            subjects_dir = self._gen_subjects_dir()
        
        outputs = self.output_spec().get()
        outputs['subject_id'] = self.inputs.subject_id
        outputs['subjects_dir'] = subjects_dir
        subject_path = os.path.join(subjects_dir, self.inputs.subject_id)
        label_list = [4 5 7 8 10 11 12 13 14 15 16 17 18 26 28 43 44 46 47 49 50 51 52 53
                      54 58 60 251 252 253 254 255]
        outputs['subcortical_surf'] = [os.path.join(subject_path, 'ascii', 'aseg_%d' %i)
                                       for i in  label_list]
        return outputs
