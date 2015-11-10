from nipype.interfaces.base import (TraitedSpec, BaseInterfaceInputSpec,
                                    File, BaseInterface, traits,
                                    CommandLineInputSpec, Directory, 
                                    CommandLine, isdefined)
from nipype.interfaces.freesurfer.base import FSTraitedSpec, FSCommand
from nipype.interfaces.matlab import MatlabCommand
from nipype.utils.filemanip import fname_presuffix, split_filename
import numpy as np
import os
from copy import deepcopy
from collections import Counter
from string import Template


class MRIsConvertInputSpec(FSTraitedSpec):
    """
    Uses Freesurfer's mris_convert to convert surface files to various formats
    """
    annot_file = File(exists=True, argstr="--annot %s",
    desc="input is annotation or gifti label data")

    parcstats_file = File(exists=True, argstr="--parcstats %s",
    desc="infile is name of text file containing label/val pairs")

    label_file = File(exists=True, argstr="--label %s",
    desc="infile is .label file, label is name of this label")

    scalarcurv_file = File(exists=True, argstr="-c %s",
    desc="input is scalar curv overlay file (must still specify surface)")

    functional_file = File(exists=True, argstr="-f %s",
    desc="input is functional time-series or other multi-frame data (must specify surface)")

    labelstats_outfile = File(exists=False, argstr="--labelstats %s",
    desc="outfile is name of gifti file to which label stats will be written")

    patch = traits.Bool(argstr="-p", desc="input is a patch, not a full surface")
    rescale = traits.Bool(argstr="-r", desc="rescale vertex xyz so total area is same as group average")
    normal = traits.Bool(argstr="-n", desc="output is an ascii file where vertex data")
    xyz_ascii = traits.Bool(argstr="-a", desc="Print only surface xyz to ascii file")
    vertex = traits.Bool(argstr="-v", desc="Writes out neighbors of a vertex in each row")

    scale = traits.Float(argstr="-s %.3f", desc="scale vertex xyz by scale")
    dataarray_num = traits.Int(argstr="--da_num %d", desc="if input is gifti, 'num' specifies which data array to use")

    talairachxfm_subjid = traits.String(argstr="-t %s", desc="apply talairach xfm of subject to vertex xyz")
    origname = traits.Str(argstr="-o %s", desc="read orig positions")

    in_file = File(exists=True, mandatory=True, position=-2, argstr='%s', desc='File to read/convert')
    out_file = File(argstr='./%s', position=-1, genfile=True, desc='output filename or True to generate one')
    #Not really sure why the ./ is necessary but the module fails without it

    out_datatype = traits.Enum("ico", "tri", "stl", "vtk", "gii", "mgh", "mgz", "asc", mandatory=True,
    desc="These file formats are supported:  ASCII:       .asc" \
    "ICO: .ico, .tri GEO: .geo STL: .stl VTK: .vtk GIFTI: .gii MGH surface-encoded 'volume': .mgh, .mgz")


class MRIsConvertOutputSpec(TraitedSpec):
    """
    Uses Freesurfer's mris_convert to convert surface files to various formats
    """
    converted = File(exists=True, desc='converted output surface')


class MRIsConvert(FSCommand):
    """
    Uses Freesurfer's mris_convert to convert surface files to various formats
    This version is a corrected version from the nipype library
    Example
    -------
    >>> import nipype.interfaces.freesurfer as fs
    >>> mris = fs.MRIsConvert()
    >>> mris.inputs.in_file = 'lh.pial'
    >>> mris.inputs.out_datatype = 'gii'
    >>> mris.run() # doctest: +SKIP
    """
    _cmd = 'mris_convert'
    input_spec = MRIsConvertInputSpec
    output_spec = MRIsConvertOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs["converted"] = os.path.abspath(self._gen_outfilename())
        return outputs

    def _gen_filename(self, name):
        if name is 'out_file':
            return self._gen_outfilename()
        else:
            return None

    def _gen_outfilename(self):
        if isdefined(self.inputs.annot_file):
            _, name, ext = split_filename(self.inputs.annot_file)
        elif isdefined(self.inputs.parcstats_file):
            _, name, ext = split_filename(self.inputs.parcstats_file)
        elif isdefined(self.inputs.label_file):
            _, name, ext = split_filename(self.inputs.label_file)
        elif isdefined(self.inputs.scalarcurv_file):
            _, name, ext = split_filename(self.inputs.scalarcurv_file)
        elif isdefined(self.inputs.functional_file):
            _, name, ext = split_filename(self.inputs.functional_file)
        elif isdefined(self.inputs.in_file):
            _, name, ext = split_filename(self.inputs.in_file)

        return name + ext + "_converted." + self.inputs.out_datatype


class Fs2TxtInputSpec(BaseInterfaceInputSpec):
    surface = File(exists=True, mandatory=True)
    out_file_triangles = File(genfile=True,
                        desc='output filename for triangles or True to generate one')
    out_file_vertices = File(genfile=True,
                        desc='output filename for vertices or True to generate one')


class Fs2TxtOutputSpec(TraitedSpec):
    vertices = File(exists=True)
    triangles = File(exists=True)


class Fs2Txt(BaseInterface):
    """
    Extract the vertice and triangle from the FreeSurfer ascii file
    """
    input_spec = Fs2TxtInputSpec
    output_spec = Fs2TxtOutputSpec

    def _get_outfilename(self):
        out_file_vertices = self.inputs.out_file_vertices
        # TODO: check if it is the right way (newpath), or use split_filename
        if not isdefined(out_file_vertices):
            out_file_vertices = fname_presuffix(self.inputs.surface,
                                      newpath=os.getcwd(),
                                      suffix='_vertices.txt',
                                      use_ext=False)
        out_file_triangles = self.inputs.out_file_triangles
        if not isdefined(out_file_triangles):
            out_file_triangles = fname_presuffix(self.inputs.surface,
                                      newpath=os.getcwd(),
                                      suffix='_triangles.txt',
                                      use_ext=False)
        return (os.path.abspath(out_file_triangles), os.path.abspath(out_file_vertices))

    def _run_interface(self, runtime):
        (out_file_triangles, out_file_vertices) = self._get_outfilename()
        with open(self.inputs.surface, 'r') as f:
            f.readline()
            nb_vert = f.readline().split(' ')[0]
            read_data = [[np.double(line.rstrip('\n').split()[0]),
                         np.double(line.rstrip('\n').split()[1]),
                         np.double(line.rstrip('\n').split()[2])] for line in f]
        a = np.array(read_data)
        vert_high = a[0:int(nb_vert), 0:3]
        tri_high = a[int(nb_vert):, 0:3]
        np.savetxt(out_file_vertices, vert_high, fmt='%.6f %.6f %.6f')
        tri_high = a[int(nb_vert):, 0:3]
        np.savetxt(out_file_triangles, tri_high, fmt='%d %d %d')
        return runtime 

    def _list_outputs(self):
        outputs = self._outputs().get()
        (out_file_triangles, out_file_vertices) = self._get_outfilename()
        outputs['vertices'] = out_file_vertices
        outputs['triangles'] = out_file_triangles
        return outputs


class Txt2OffInputSpec(BaseInterfaceInputSpec):
    vertices = File(exists=True, mandatory=True)
    triangles = File(exists=True, mandatory=True)
    out_file_off = File(genfile=True,
                        desc='output filename for off file or True to generate one')
    

class Txt2OffOutputSpec(TraitedSpec):
    surface_off = File(exists=True)


class Txt2Off(BaseInterface):
    """
    Convert the .txt files into .off files for use in remesher 
    """
    input_spec = Txt2OffInputSpec
    output_spec = Txt2OffOutputSpec

    def _get_outfilename(self):
        out_file_off = self.inputs.out_file_off
        if not isdefined(out_file_off):
            out_file_off= fname_presuffix(self.inputs.vertices,
                                      newpath=os.getcwd(),
                                      suffix='.off',
                                      use_ext=False)
        return os.path.abspath(out_file_off)

    def _run_interface(self, runtime):
        vert = np.loadtxt(self.inputs.vertices)
        tri = np.loadtxt(self.inputs.triangles)
        out_file_off = self._get_outfilename()
        with open(out_file_off, 'w') as f: 
            f.write('OFF\n') 
            f.write('{} {} {}\n'.format(int(vert.shape[0]), int(tri.shape[0]), 0)) 
        with open(out_file_off, 'a') as f:
            np.savetxt(f, vert, fmt='%.6f')
            np.savetxt(f, np.hstack([np.ones((tri.shape[0],1))*3, tri]), fmt='%d')
        return runtime 

    def _list_outputs(self):
        outputs = self._outputs().get()
        out_file_off = self._get_outfilename()
        outputs['surface_off'] = out_file_off
        return outputs


class RemesherInputSpec(CommandLineInputSpec):
    in_file = File(desc="Input surface", argstr='%s', exists=True, 
                   mandatory=True, position=0)
    out_file = File(desc="Remeshed surface", argstr='%s',
                    genfile=True, position=1)
    scripts_directory = traits.Str(desc="scripts directory", exists=True)


class RemesherOutputSpec(TraitedSpec):
    out_file = File(desc = "Remeshed surface", exists = True)


class Remesher(CommandLine):
    """
    Wrapper for remesher command 
    For input rh_high.txt will return output rh_low.txt
    """
    input_spec = RemesherInputSpec
    output_spec = RemesherOutputSpec
    _cmd = os.getcwd()

    def __init__(self, *args, **kwargs):
        super(Remesher, self).__init__(*args, **kwargs)
        self._cmd = self.get_scripts_dir() + "/remesher-mac/cmdremesher/cmdremesher"
    
    def get_scripts_dir(self):
        if isdefined(self.inputs.scripts_directory):
            return self.inputs.scripts_directory
        else:
            return os.getcwd()

    def _gen_filename(self, name):
        if name is 'out_file':
            return self._gen_outfilename()
        else:
            return None
            
    def _gen_outfilename(self):
        if isdefined(self.inputs.out_file):
            return self.inputs.out_file
        else:
            _, name, ext = split_filename(self.inputs.in_file)
            return name[:2] + "_decimated" + ext 

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_file'] = os.path.abspath(self._gen_outfilename())
        return outputs


class Off2TxtInputSpec(BaseInterfaceInputSpec):
    surface_off = File(exists=True, mandatory=True)
    out_file_vertices_txt = File(genfile=True,
                        desc='output filename for off file or True to generate one')
    out_file_triangles_txt = File(genfile=True,
                        desc='output filename for off file or True to generate one')
    

class Off2TxtOutputSpec(TraitedSpec):
    vertices_txt = File(exists=True)
    triangles_txt = File(exists=True)


class Off2Txt(BaseInterface):
    """
    Convert the .off files into .txt files for use in remesher 
    """
    input_spec = Off2TxtInputSpec
    output_spec = Off2TxtOutputSpec

    def _get_outfilename(self):
        out_file_vertices_txt = self.inputs.out_file_vertices_txt
        out_file_triangles_txt = self.inputs.out_file_triangles_txt
        if not isdefined(out_file_triangles_txt):
            out_file_triangles_txt = fname_presuffix(self.inputs.surface_off,
                                      newpath=os.getcwd(),
                                      suffix='_triangles.txt',
                                      use_ext=False)
        if not isdefined(out_file_vertices_txt):
            out_file_vertices_txt = fname_presuffix(self.inputs.surface_off,
                                      newpath=os.getcwd(),
                                      suffix='_vertices.txt',
                                      use_ext=False)
        return (os.path.abspath(out_file_triangles_txt), os.path.abspath(out_file_vertices_txt))

    def _run_interface(self, runtime):
        out_file_triangles_txt, out_file_vertices_txt = self._get_outfilename()
        with open(self.inputs.surface_off) as f:
            f.readline()
            num = f.readline().split(' ')
            vert = np.loadtxt(self.inputs.surface_off, skiprows=2, usecols=(0,1,2))  
            vert = vert[:int(num[0]), :]
            tri = np.loadtxt(self.inputs.surface_off, skiprows=int(num[0])+2, usecols=(1,2,3))  
        np.savetxt(out_file_vertices_txt, vert, fmt='%.4f')
        np.savetxt(out_file_triangles_txt, tri, fmt='%d')
        return runtime 

    def _list_outputs(self):
        outputs = self._outputs().get()
        out_file_triangles_txt, out_file_vertices_txt = self._get_outfilename()        
        outputs['vertices_txt'] = out_file_vertices_txt
        outputs['triangles_txt'] = out_file_triangles_txt
        return outputs


class RegionMappingInputSpect(BaseInterfaceInputSpec):
    aparc_annot = File(exists = True, mandatory = True)
    ref_table = File(exists = True, mandatory = True)
    vertices_downsampled = File(exists = True, mandatory = True)
    triangles_downsampled = File(exists = True, mandatory = True)
    vertices = File(exists = True, mandatory = True)
    scripts_directory = Directory(desc="scripts directory", exists=True)
    out_file = File(desc="Region mapping")

class RegionMappingOutputSpect(TraitedSpec):
    out_file = File(exists = True)


class RegionMapping(BaseInterface):
    """
    Generate a first region_mapping txt file using the region_mapping_2 matlab function  
    """
    input_spec = RegionMappingInputSpect
    output_spec = RegionMappingOutputSpect

    def _gen_outfilename(self):
        if isdefined(self.inputs.out_file):
            return self.inputs.out_file
        else:
            _, name, ext = split_filename(self.inputs.vertices)
            return os.path.abspath(name + "_region_mapping" + ext)
    
    def get_scripts_dir(self):
        if isdefined(self.inputs.scripts_directory):
            return self.inputs.scripts_directory
        else:
            return os.getcwd()

    def _run_interface(self, runtime):
        d = dict(vertices_downsampled = self.inputs.vertices_downsampled,
                 triangles_downsampled = self.inputs.triangles_downsampled,
                 vertices = self.inputs.vertices,
                 ref_tables = self.inputs.ref_table,
                 aparc_annot = self.inputs.aparc_annot,
                 out_file = self._gen_outfilename(),
                 scripts_dir = os.path.abspath(self.get_scripts_dir())
                 )
        script = Template("""
            vertices_downsampled = '$vertices_downsampled';
            triangles_downsampled = '$triangles_downsampled';
            vertices = '$vertices';
            ref_tables = '$ref_tables';
            aparc_annot = '$aparc_annot';
            out_file = '$out_file';
            addpath($scripts_dir);
            $scripts_dir/region_mapping(vertices_downsampled, triangles_downsampled,
                             vertices, ref_tables, aparc_annot, out_file); 
            quit;
            """).safe_substitute(d)
        mlab = MatlabCommand(script=script, mfile=True)
        result = mlab.run()
        return result.runtime

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = self._gen_outfilename()
        return outputs



class CorrectRegionMappingInputSpec(BaseInterfaceInputSpec):
    vertices = File(exists=True, mandatory=True)
    triangles = File(exists=True, mandatory=True)
    texture = File(exists=True, mandatory=True)
    region_mapping_corr = traits.Float(exists=True, mandatory=True)
    out_file_region_mapping_corr = File(genfile=True,
                                        desc='output filename for off file or True to generate one')


class CorrectRegionMappingOutputSpec(TraitedSpec):
    texture_corrected= File(exists=True)


class CorrectRegionMapping(BaseInterface):
    """
    Correct the region mapping
    """
    input_spec = CorrectRegionMappingInputSpec
    output_spec = CorrectRegionMappingOutputSpec

    def _get_outfilename(self):
        out_file_region_mapping_corr = self.inputs.out_file_region_mapping_corr
        if not isdefined(out_file_region_mapping_corr):
            out_file_region_mapping_corr= fname_presuffix(self.inputs.in_file,
                                      newpath=os.getcwd(),
                                      suffix='_corr',
                                      use_ext=True)
        return os.path.abspath(out_file_region_mapping_corr)

    def _run_interface(self, runtime):
        texture = loadtxt(self.inputs.texture)
        vert = loadtxt(self.inputs.vertices)
        trian = loadtxt(self.inputs.triangles)
        out_file_region_mapping_corr = self._get_outfilename()
        for _ in range(10):
           new_texture = deepcopy(texture)
           labels = np.unique(texture)
           for ilab in labels:
               iverts = (np.nonzero(texture==ilab)[0]).tolist()
               if len(iverts)>0:
                   for inode in iverts:
                       iall = trian[np.nonzero(trian==inode)[0]].flatten().tolist()
                       ineig = np.unique(filter(lambda x: x!=inode, iall)).astype('int')
                       ivals = np.array(Counter(texture[ineig]).most_common()).astype('int') 
                       if ivals[np.nonzero(ivals[:,0]==ilab), 1].shape[1]==0:
                           new_texture[inode] = Counter(texture[ineig]).most_common(1)[0][0]
                       elif ivals[np.nonzero(ivals[:,0] == ilab), 1][0,0] < region_mapping_corr * len(ineig):
                           new_texture[inode] = Counter(texture[ineig]).most_common(1)[0][0]
           texture = deepcopy(new_texture) 
        
        np.savetxt(out_file_region_mapping_corr, new_texture)
        return runtime 

    def _list_outputs(self):
        outputs = self._outputs().get()
        out_file_region_mapping_corr = self._get_outfilename()
        outputs['texture_corrected'] = out_file_region_mapping_corr
        return outputs


class CheckRegionMappingInputSpect(CommandLineInputSpec):
    vertices = File(argstr ='%s', exists = True, 
                        mandatory = True, position = 0)
    triangles = File(argstr ='%s', exists = True, 
                         mandatory = True, position = 1)
    region_mapping = File(argstr ='%s', exists = True, 
                              mandatory = True, position = 2)


class CheckRegionMappingOutputSpect(TraitedSpec):
    region_mapping = File(exists=True)


class CheckRegionMapping(CommandLine):
    input_spec = CheckRegionMappingInputSpect
    output_spec = CheckRegionMappingOutputSpect
    _cmd = 'python /home/tim/Work/Models/processing/scripts_nathan/check_region_mapping_2.py'
    #_cmd = scripts_dir + "/check_region_mapping_2.py"
    _terminal_output = 'stream'
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['region_mapping_low'] = self.inputs.region_mapping_low
        return outputs


class ReunifyBothHemisphereInputSpec(BaseInterfaceInputSpec):
    vertices = traits.List(traits=File, exists=True, mandatory=True, minlen=2, maxlen=2, 
                    desc='right and left hemispheres for vertices')
    triangles = traits.List(traits=File, exists=True, mandatory=True, minlen=2, maxlen=2, 
                    desc='right and left hemispheres for triangles')
    textures = traits.List(traits=File, exists=True, mandatory=True, minlen=2, maxlen=2, 
                    desc='right and left hemispheres for texture')
    out_file_vertices = File(genfile=True,
                        desc='output filename for vertices file or True to generate one')
    out_file_triangles = File(genfile=True,
                        desc='output filename for triangles file or True to generate one')
    out_file_texture = File(genfile=True,
                        desc='output filename for texture file or True to generate one')


class ReunifyBothHemisphereOutputSpec(TraitedSpec):
    texture = File(exists=True)
    vertices = File(exists=True)
    triangles = File(exists=True)


class ReunifyBothHemisphere(BaseInterface):
    """
    Reunify both hemispheres (triangles, vertices and region mapping) 
    """
    input_spec = ReunifyBothHemisphereInputSpec
    output_spec = ReunifyBothHemisphereOutputSpec

    def _get_outfilename(self):
        out_file_vertices = self.inputs.out_file_vertices
        out_file_triangles = self.inputs.out_file_triangles
        out_file_texture = self.inputs.out_file_texture
        if not isdefined(out_file_triangles):
            out_file_triangles= fname_presuffix(self.inputs.triangles_rh,
                                      newpath=os.getcwd(),
                                      suffix='_triangles',
                                      use_ext=True)
        if not isdefined(out_file_vertices):
            out_file_vertices= fname_presuffix(self.inputs.vertices_rh,
                                      newpath=os.getcwd(),
                                      suffix='_vertices',
                                      use_ext=True)
        if not isdefined(out_file_texture):
            out_file_texture= fname_presuffix(self.inputs.texture_rh,
                                      newpath=os.getcwd(),
                                      suffix='_texture',
                                      use_ext=True)
        return (os.path.abspath(out_file_triangles), 
                os.path.abspath(out_file_vertices),
                os.path.abspath(out_file_texture))

    def _run_interface(self, runtime):
        (out_file_triangles, out_file_vertices, out_file_texture) = self._get_outfilename()
        lh_reg_map = np.loadtxt(self.inputs.texture[0])
        rh_reg_map = np.loadtxt(self.inputs.texture[1])
        lh_vert = np.loadtxt(self.inputs.vertices[0])
        rh_vert = np.loadtxt(self.inputs.vertices[1])
        lh_trian = np.loadtxt(self.inputs.triangles[0])
        rh_trian = np.loadtxt(self.inputs.triangels[1])
        vertices = vstack([lh_vert, rh_vert])
        triangles = vstack([lh_trian,  rh_trian + lh_vert.shape[0]])
        region_mapping = hstack([lh_reg_map, rh_reg_map])
        np.savetxt(out_file_texture, region_mapping, fmt='%d', newline=" ")
        np.savetxt(out_file_vertices, vertices, fmt='%.2f')
        np.savetxt(out_file_triangles, triangles, fmt='%d %d %d')
        return runtime 

    def _list_outputs(self):
        outputs = self._outputs().get()
        (out_file_triangles, out_file_vertices, out_file_texture) = self._get_outfilename()
        outputs['texture'] = out_file_texture
        outputs['vertices'] = out_file_vertices
        outputs['triangles'] = out_file_triangles
        return outputs


class Aseg2SrfInputSpec(CommandLineInputSpec):
    subject_id = traits.Str(desc="FreeSurfer Subject Id", argstr='-s %s',
                      exists=True, mandatory=True)
    subjects_dir = Directory(desc='path to subject directory', argstr='-sd %s',
                             exists=True, mandatory=True)
    # list in Traits
    label_list = traits.List([4, 5, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 26, 28, 43,
                      44, 46, 47, 49, 50, 51, 52, 53, 54, 58, 60, 251, 252, 253,
                      254, 255])


class Aseg2SrfOutputSpec(TraitedSpec):
    subcortical_surf_list = traits.List(File(desc="Output subcortical surfaces",
                                 exists=True))


class Aseg2Srf(CommandLine):
    input_spec = Aseg2SrfInputSpec
    output_spec = Aseg2SrfOutputSpec
    # TODO: how to use a command which is not in the path?
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
        outputs['subcortical_surf'] = [os.path.join(subject_path, 'ascii', 'aseg_%d' %i)
                                       for i in  self.label_list]
        return outputs









class ListSubcorticalInputSpect(CommandLineInputSpec):
    in_file = File(argstr ='%s', 
                        exists = True, 
                        mandatory = True,
                        position = 0)

class ListSubcorticalOutputSpect(TraitedSpec):
    triangles = File(exists=True)
    vertices = File(exists=True)

class ListSubcortical(CommandLine):
    input_spec = ListSubcorticalInputSpect
    output_spec = ListSubcorticalOutputSpect
    _cmd = 'python /home/tim/Work/Models/processing/scripts_nathan/list_subcortical_2.py'
    #_cmd = scripts_dir + "/list_subcortical_2.py"
    _terminal_output = 'stream'
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['triangles'] = os.path.abspath(fname_presuffix("",  prefix=self.inputs.in_file[-12:-4], suffix='_tri.txt'))
        outputs['vertices'] = os.path.abspath(fname_presuffix("",  prefix=self.inputs.in_file[-12:-4], suffix='_vert.txt'))
        return outputs





