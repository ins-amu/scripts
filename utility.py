from nipype.interfaces.base import TraitedSpec, BaseInterfaceInputSpec, File, BaseInterface
import numpy as np
import os
from copy import deepcopy
from collections import Counter

class ExtractHighInputSpec(BaseInterfaceInputSpec):
    surface = File(exists=True, mandatory=True)


class ExtractHighOutputSpec(TraitedSpec):
    vertices_high = File(exists=True)
    triangles_high = File(exists=True)


class ExtractHigh(BaseInterface):
    """
    Extract the vertice and triangle from the FreeSurfer ascii file
    """
    input_spec = ExtractHighInputSpec
    output_spec = ExtractHighOutputSpec

    def _run_interface(self, runtime):
        with open(self.inputs.surface, 'r') as f:
            f.readline()
            nb_vert = f.readline().split(' ')[0]
            read_data = [[np.double(line.rstrip('\n').split()[0]),
                         np.double(line.rstrip('\n').split()[1]),
                         np.double(line.rstrip('\n').split()[2])] for line in f]
        a = np.array(read_data)
        vert_high = a[0:int(nb_vert), 0:3]
        tri_high = a[int(nb_vert):, 0:3]
        np.savetxt('vertices_high.txt', vert_high, fmt='%.6f %.6f %.6f')
        tri_high = a[int(nb_vert):, 0:3]
        np.savetxt('triangles_high.txt', tri_high, fmt='%d %d %d')
        return runtime 

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['vertices_high'] = os.path.abspath('vertices_high.txt') 
        return outputs


class Txt2OffInputSpec(BaseInterfaceInputSpec):
    vertices = File(exists=True, mandatory=True)
    triangles = File(exists=True, mandatory=True)


class Txt2OffOutputSpec(TraitedSpec):
    surface_high_off = File(exists=True)


class Txt2Off(BaseInterface):
    """
    Convert the .txt files into .off files for use in remesher 
    """
    input_spec = Txt2OffInputSpec
    output_spec = Txt2OffHighOutputSpec

    def _run_interface(self, runtime):
        vert = np.loadtxt(self.inputs.vertices)
        tri = np.loadtxt(self.inputs.triangles)
        with open('surface_high.off', 'w') as f: f.write('OFF\n') 
            f.write('{} {} {}\n'.format(int(vert.shape[0]), int(tri.shape[0]), 0)) 
        with open('surface_high.off', 'a') as f:
            np.savetxt(f, vert, fmt='%.6f')
            np.savetxt(f, np.hstack([np.ones((tri.shape[0],1))*3, tri]), fmt='%d')
        return runtime 

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['surface_high_off'] = os.path.abspath('surface_high.off') 
        return outputs


class Off2TxtInputSpec(BaseInterfaceInputSpec):
    surface_low_off = File(exists=True, mandatory=True)


class Off2TxtOutputSpec(TraitedSpec):
    vertices = File(exists=True)
    triangles = File(exists=True)


class Off2Txt(BaseInterface):
    """
    Convert the .off files into .txt files for use in remesher 
    """
    input_spec = Off2TxtInputSpec
    output_spec = Off2TxtHighOutputSpec

    def _run_interface(self, runtime):
        with open(self.inputs.surface) as f:
            f.readline()
            num = f.readline().split(' ')
            vert = np.loadtxt(surface, skiprows=2, usecols=(0,1,2))  
            vert = vert[:int(num[0]), :]
            tri = np.loadtxt(surface, skiprows=int(num[0])+2, usecols=(1,2,3))  
        np.savetxt('vertices_low.off', vert, fmt='%.4f')
        np.savetxt('triangles_low.off', tri, fmt='%d')
        return runtime 

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['vertices_low_off'] = os.path.abspath('vertices_low.off') 
        outputs['triangles_low_off'] = os.path.abspath('triangles_low.off') 
        return outputs


class CorrectRegionMappingInputSpec(BaseInterfaceInputSpec):
    vertices = File(exists=True, mandatory=True)
    triangles = File(exists=True, mandatory=True)
    texture = File(exists=True, mandatory=True)
    region_mapping_corr = Float(exists=True, mandatory=True)


class CorrectRegionMappingOutputSpec(TraitedSpec):
    texture_corrected= File(exists=True)


class CorrectRegionMapping(BaseInterface):
    """
    Correct the region mapping
    """
    input_spec = CorrectRegionMappingInputSpec
    output_spec = CorrectRegionMappingOutputSpec

    def _run_interface(self, runtime):
        texture = loadtxt(self.inputs.texture)
        vert = loadtxt(self.inputs.vertices)
        trian = loadtxt(self.inputs.triangles)
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
        
        np.savetxt('region_mapping_low.txt', new_texture)
        return runtime 

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['texture_corrected'] = os.path.abspath('region_mapping_low.txt') 
        return outputs

class ReunifyBothHemisphereInputSpec(BaseInterfaceInputSpec):
    vertices_rh = File(exists=True, mandatory=True)
    vertices_lh = File(exists=True, mandatory=True)
    triangles_rh = File(exists=True, mandatory=True)
    triangles_lh = File(exists=True, mandatory=True)
    texture_rh = File(exists=True, mandatory=True)
    texture_lh = File(exists=True, mandatory=True)


class ReunifyBothHemisphereOutputSpec(TraitedSpec):
    texture = File(exists=True)
    vertices = File(exists=True)
    triangles = File(exists=True)


class ReunifyBothHemisphere(BaseInterface):
    """
    Reunify both hemispheres (triangles, vertices and region mapping) 
    """
    input_spec = CorrectRegionMappingInputSpec
    output_spec = CorrectRegionMappingOutputSpec

    def _run_interface(self, runtime):
        lh_reg_map = np.loadtxt(self.inputs.texture_lh)
        lh_vert = np.loadtxt(self.inputs.vertices_lh)
        lh_trian = np.loadtxt(self.inputs.triangles_lh)
        rh_reg_map = np.loadtxt(self.inputs.texture_rh)
        rh_vert = np.loadtxt(self.inputs.texture_rh)
        rh_trian = np.loadtxt(self.inputs.texture_rh)
        vertices = vstack([lh_vert, rh_vert])
        triangles = vstack([lh_trian,  rh_trian + lh_vert.shape[0]])
        region_mapping = hstack([lh_reg_map, rh_reg_map])
        np.savetxt('region_mapping.txt', region_mapping, fmt='%d', newline=" ")
        np.savetxt('vertices.txt', vertices, fmt='%.2f')
        np.savetxt('triangles.txt', triangles, fmt='%d %d %d')
        return runtime 

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['texture'] = os.path.abspath('region_mapping.txt') 
        outputs['vertices'] = os.path.abspath('vertices.txt') 
        outputs['triangles'] = os.path.abspath('triangles.txt') 
        return outputs


class Aseg2SrfInputSpec(CommandLineSpec):
    subject_id = File(desc = "Subject FreeSurfer Id", argstr = '%d',
            exists = True, mandatory = True)


class Aseg2SrfOutputSpec(TraitedSpec):
    subcortical_surf_list = File(desc = "Output subcortical surfaces", exists = True)


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

class RemesherInputSpec(CommandLineInputSpec):
    in_file = File(desc = "Input surface", 
                   argstr ='%s', 
                   exists = True, 
                   mandatory = True,
                   position = 0)
    out_file = File(desc = "Remeshed surface",
                    argstr ='%s',
                    exists = True,
                    genfile=True,
                    position = 1)

class RemesherOutputSpec(TraitedSpec):
    out_file = File(desc = "Remeshed surface", exists = True)

class Remesher(CommandLine):
    """
    Wrapper for remesher command 
    For input rh_high.txt will return output rh_low.txt
    """
    input_spec = RemesherInputSpec
    output_spec = RemesherOutputSpec
    _cmd = "/home/tim/Work/Models/processing/scripts_nathan/remesher/cmdremesher/cmdremesher"
    #_cmd = scripts_dir + "/remesher/cmdremesher/cmdremesher"
    def _gen_filename(self, name):
        if name is 'out_file':
            return self._gen_outfilename()
        else:
            return None

    def _gen_outfilename(self):
        if isdefined(self.inputs.in_file):
            _, name, ext = split_filename(self.inputs.in_file)
            return name[:2] + "_low" + ext 

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_file'] = os.path.abspath(self._gen_outfilename())
        return outputs


class RegionMappingInputSpect(BaseInterfaceInputSpec):
    rl = traits.Str(mandatory=True,desc=("right or left hemisphere"))
    aparc_annot = File(exists = True, mandatory = True)
    ref_tables = File(exists = True, mandatory = True)
    vertices_low = File(exists = True, mandatory = True)
    triangles_low = File(exists = True, mandatory = True)
    vertices_high = File(exists = True, mandatory = True)

class RegionMappingOutputSpect(TraitedSpec):
    out_file = File(exists = True)

class RegionMapping(BaseInterface):
    """
    Generate a first region_mapping txt file using the region_mapping_2 matlab function  
    """
    input_spec = RegionMappingInputSpect
    output_spec = RegionMappingOutputSpect

    def _run_interface(self, runtime):
        d = dict(rl = self.inputs.rl,
                 vertices_low = self.inputs.vertices_low,
                 triangles_low = self.inputs.triangles_low,
                 vertices_high = self.inputs.vertices_high,
                 ref_tables = self.inputs.ref_tables,
                 aparc_annot = self.inputs.aparc_annot
                 # scripts_dir = scripts_dir
                 )
        script = Template("""
            rl = '$rl'
            vertices_low = '$vertices_low';
            triangles_low = '$triangles_low';
            vertices_high = '$vertices_high';
            ref_tables = '$ref_tables';
            aparc_annot = '$aparc_annot';
            addpath('/home/tim/Work/Models/processing/scripts_nathan');
            region_mapping_2(rl, vertices_low, triangles_low, vertices_high, ref_tables, aparc_annot); 
            quit;
            """).safe_substitute(d)
        ##changes needed for addpath (specific architecture) ->addpath('$scipts_dir');
        mlab = MatlabCommand(script=script, mfile=True)
        result = mlab.run()
        return result.runtime

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = os.path.abspath(fname_presuffix("",  prefix=self.inputs.rl, suffix='_region_mapping_low_not_corrected.txt'))
        return outputs


class CheckRegionMappingInputSpect(CommandLineInputSpec):
    vertices_low = File(argstr ='%s', 
                        exists = True, 
                        mandatory = True,
                        position = 0)
    triangles_low = File(argstr ='%s', 
                         exists = True, 
                         mandatory = True,
                         position = 1)
    region_mapping_low = File(argstr ='%s', 
                              exists = True, 
                              mandatory = True,
                              position = 2)

class CheckRegionMappingOutputSpect(TraitedSpec):
    region_mapping_low = File(exists=True)

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
    origname = traits.String(argstr="-o %s", desc="read orig positions")

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


