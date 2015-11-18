from nipype.interfaces.base import (TraitedSpec, BaseInterfaceInputSpec,
                                    File, BaseInterface, traits,
                                    CommandLineInputSpec, Directory,
                                    CommandLine, isdefined)
from nipype.interfaces.freesurfer.base import FSTraitedSpec, FSCommand
from nipype.interfaces.matlab import MatlabCommand
from nipype.utils.filemanip import fname_presuffix, split_filename
import numpy as np
import os, sys
from copy import deepcopy
from collections import Counter
from string import Template


class GrabDataInputSpec(BaseInterfaceInputSpec):
    subject_id = traits.Str(exists=True, mandatory=True)
    freesurfer_directory = Directory(exists=True, mandatory=True)
    scripts_directory = Directory(exists=True)
    hemi = traits.Str(exists=True, mandatory=True)
    out_file_pial = File(desc='output filename for pial file')
    out_file_annot = File(desc='output filename for annot file')
    out_file_ref_table = File(desc='output filename for ref table file')


class GrabDataOutputSpec(TraitedSpec):
    pial = File(desc='output filename for pial file')
    annot = File(desc='output filename for annot file')
    ref_table = File(desc='output filename for ref table file')


class GrabData(BaseInterface):
    """
    grab the files produced by freesurfer.
    """
    input_spec = GrabDataInputSpec
    output_spec = GrabDataOutputSpec

    def _get_scripts_directory(self):
        if isdefined(self.inputs.scripts_directory):
            out_file_scripts = self.inputs.scripts_directory
        else:
            out_file_scripts = os.getcwd()
        return os.path.abspath(out_file_scripts)

    def _run_interface(self, runtime):
        return runtime

    def _list_outputs(self):
        path_fs_subj = os.path.join(os.path.abspath(self.inputs.freesurfer_directory), self.inputs.subject_id)
        out_file_pial = os.path.join(path_fs_subj, 'surf', self.inputs.hemi + '.pial')
        out_file_annot = os.path.join(path_fs_subj, 'label', self.inputs.hemi + '.aparc.annot')
        out_file_ref_table = os.path.join(self._get_scripts_directory(), self.inputs.hemi + '_ref_table.txt')
        outputs = self._outputs().get()
        outputs['pial'] = out_file_pial
        outputs['annot'] = out_file_annot
        outputs['ref_table'] = out_file_ref_table
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
    origname = traits.Str(argstr="-o %s", desc="read orig positions")

    in_file = File(exists=True, mandatory=True, position=-2, argstr='%s', desc='File to read/convert')
    out_file = File(argstr='./%s', position=-1, genfile=True, desc='output filename or True to generate one')
    # Not really sure why the ./ is necessary but the module fails without it

    out_datatype = traits.Enum("ico", "tri", "stl", "vtk", "gii", "mgh", "mgz", "asc", mandatory=True,
                               desc="These file formats are supported:  ASCII: .asc"
                                    "ICO: .ico, .tri GEO: .geo STL: .stl VTK: .vtk GIFTI: .gii"
                                    "MGH surface-encoded 'volume': .mgh, .mgz")


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
            out_file_off = fname_presuffix(self.inputs.vertices,
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
            np.savetxt(f, np.hstack([np.ones((tri.shape[0], 1)) * 3, tri]), fmt='%d')
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
    out_file = File(desc="Remeshed surface", exists=True)


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
        if sys.platform=='darwin':
            self._cmd = self.get_scripts_dir() + "/remesher-mac/cmdremesher/cmdremesher"
        elif sys.platform=='linux2':
            self._cmd = self.get_scripts_dir() + "/remesher-mac/cmdremesher/cmdremesher"
        else:
            raise NotImplementedError, 'only working on mac or linux'

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
            vert = np.loadtxt(self.inputs.surface_off, skiprows=2, usecols=(0, 1, 2))
            vert = vert[:int(num[0]), :]
            tri = np.loadtxt(self.inputs.surface_off, skiprows=int(num[0]) + 2, usecols=(1, 2, 3))
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
    aparc_annot = File(exists=True, mandatory=True)
    ref_table = File(exists=True, mandatory=True)
    vertices_downsampled = File(exists=True, mandatory=True)
    triangles_downsampled = File(exists=True, mandatory=True)
    vertices = File(exists=True, mandatory=True)
    # try with getcwd
    scripts_directory = Directory(desc="scripts directory", exists=True)
    out_file = File(desc="Region mapping")


class RegionMappingOutputSpect(TraitedSpec):
    out_file = File(exists=True)


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
        d = dict(vertices_downsampled=self.inputs.vertices_downsampled,
                 triangles_downsampled=self.inputs.triangles_downsampled,
                 vertices=self.inputs.vertices,
                 ref_tables=self.inputs.ref_table,
                 aparc_annot=self.inputs.aparc_annot,
                 out_file=self._gen_outfilename(),
                 scripts_dir=os.path.abspath(self.get_scripts_dir())
                 )
        script = Template("""
            vertices_downsampled = '$vertices_downsampled';
            triangles_downsampled = '$triangles_downsampled';
            vertices = '$vertices';
            ref_tables = '$ref_tables';
            aparc_annot = '$aparc_annot';
            out_file = '$out_file';
            addpath('$scripts_dir');
            region_mapping(vertices_downsampled, triangles_downsampled, vertices, ref_tables, aparc_annot, out_file); 
            exit;
            """).substitute(d)
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
    region_mapping_corr = traits.Float(usedefault=True)
    out_file_region_mapping_corr = File(desc='output filename for corrected region mapping')


class CorrectRegionMappingOutputSpec(TraitedSpec):
    texture_corrected = File(exists=True)


class CorrectRegionMapping(BaseInterface):
    """
    Correct the region mapping
    """
    input_spec = CorrectRegionMappingInputSpec
    output_spec = CorrectRegionMappingOutputSpec
    region_mapping_corr = traits.Float(0.42)

    def _get_outfilename(self):
        out_file_region_mapping_corr = self.inputs.out_file_region_mapping_corr
        if not isdefined(out_file_region_mapping_corr):
            out_file_region_mapping_corr = fname_presuffix(self.inputs.texture, newpath=os.getcwd(),
                                                           suffix='_corr', use_ext=True)
        return os.path.abspath(out_file_region_mapping_corr)

    def _run_interface(self, runtime):
        texture = np.loadtxt(self.inputs.texture)
        # vert = np.loadtxt(self.inputs.vertices)
        trian = np.loadtxt(self.inputs.triangles)
        out_file_region_mapping_corr = self._get_outfilename()
        for _ in range(10):
            new_texture = deepcopy(texture)
            labels = np.unique(texture)
            for ilab in labels:
                iverts = (np.nonzero(texture == ilab)[0]).tolist()
                if len(iverts) > 0:
                    for inode in iverts:
                        iall = trian[np.nonzero(trian == inode)[0]].flatten().tolist()
                        ineig = np.unique(filter(lambda x: x != inode, iall)).astype('int')
                        ivals = np.array(Counter(texture[ineig]).most_common()).astype('int')
                        if ivals[np.nonzero(ivals[:, 0] == ilab), 1].shape[1] == 0:
                            new_texture[inode] = Counter(texture[ineig]).most_common(1)[0][0]
                        # TODO replace by inputs region_mapping_corr
                        elif ivals[np.nonzero(ivals[:, 0] == ilab), 1][0, 0] < 0.42 * len(ineig):
                            new_texture[inode] = Counter(texture[ineig]).most_common(1)[0][0]
            texture = deepcopy(new_texture)

        np.savetxt(out_file_region_mapping_corr, new_texture)
        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        out_file_region_mapping_corr = self._get_outfilename()
        outputs['texture_corrected'] = out_file_region_mapping_corr
        return outputs


class CheckRegionMappingInputSpect(BaseInterfaceInputSpec):
    vertices = File(exists=True, mandatory=True)
    triangles = File(exists=True, mandatory=True)
    region_mapping = File(exists=True, mandatory=True)
    scripts_directory = Directory(desc="scripts directory", exists=True)
    check = traits.Bool(mandatory=True, desc="check the region mapping")
    display = traits.Bool(mandatory=True, desc="display the region mapping")
    out_file = File(desc='out file name')


class CheckRegionMappingOutputSpect(TraitedSpec):
    region_mapping = File(exists=True)


class CheckRegionMapping(BaseInterface):
    input_spec = CheckRegionMappingInputSpect
    output_spec = CheckRegionMappingOutputSpect

    # TODO: more generic
    # _terminal_output = 'stream'

    def _get_outfilename(self):
        if isdefined(self.inputs.out_file):
            return self.inputs.out_file
        else:
            _, name, ext = split_filename(self.inputs.vertices)
            return os.path.abspath(name + "_checked_region_mapping" + ext)

    def get_scripts_dir(self):
        if isdefined(self.inputs.scripts_directory):
            return self.inputs.scripts_directory
        else:
            return os.getcwd()

    def breadth_first_search(self, iposi, itrian, ilab, texture):
        queue = [iposi]
        V = [iposi]
        while len(queue) > 0:
            iQ = queue.pop()
            iedges = list(itrian[np.argwhere(itrian == iQ)[:, 0]].flatten())
            while len(iedges) > 0:
                ineig = iedges.pop()
                if ineig not in V and texture[ineig] == ilab:
                    V.append(ineig)
                    queue.append(ineig)
        return V

    def calculate_connected(self, texture, vert, trian):
        """
        find if the regions are connected components using Breadth-first seach

        :param texture: np array
        :param vert: list
        :param trian: list
        """
        labels = np.unique(texture)
        res = []
        for ilab in labels:
            ipos = np.nonzero(texture == ilab)
            ivert = vert[ipos]
            itrian = []
            for itri in np.nonzero(texture == ilab)[0].tolist():
                itrian.extend(trian[np.nonzero(trian == itri)[0]])
            itrian = np.array(itrian).astype('int')
            V = self.breadth_first_search(ipos[0][0], itrian, ilab, texture)
            res.append((ilab, ivert.shape[0] - len(V)))
        return res

    def find_both_components(self, texture, vert, trian, ilab):
        """
        find the two subgraphs

        :param texture: np array
        :param vert: list
        :param trian: list
        :param ilab: list
        """
        ipos = np.nonzero(texture == ilab)
        # ivert = vert[ipos]
        itrian = []
        for itri in ipos[0].tolist():
            itrian.extend(trian[np.nonzero(trian == itri)[0]])
        itrian = np.array(itrian).astype('int')
        # first region
        V1 = self.breadth_first_search(ipos[0][0], itrian, ilab, texture)
        # second region
        istart = 1
        while ipos[0][istart] in V1:
            istart += 1
        V2 = self.breadth_first_search(ipos[0][istart], itrian, ilab, texture)
        return (V1, V2)

    def correct_sub_region(self, texture, trian, Vw):
        """
        correct the region mapping for the chosen component

        :param texture: np array
        :param trian: list
        :param Vw: list
        """
        new_texture = deepcopy(texture)
        icount = 0
        while len(Vw) > 0:
            iVw = Vw.pop()
            itrian = trian[np.nonzero(trian == iVw)[0]].flatten().astype('int').tolist()
            ir = filter(lambda x: new_texture[x] != new_texture[iVw], itrian)
            if len(ir) > 0:
                new_texture[iVw] = new_texture[Counter(ir).most_common(1)[0][0]]
            else:
                if icount < 50:
                    Vw.insert(0, iVw)
                    icount += 1
                else:
                    # TODO: good error message
                    print('error in correction')
                    import pdb
                    pdb.set_trace()
        return new_texture

    def check_region_mapping(self, texture, vert, trian, ilab):
        """
        drawing the region

        :param texture: np array
        :param vert: list
        :param trian: list
        :param ilab: list
        """
        ipos = np.nonzero(texture == ilab)
        itrian = []
        for itri in ipos[0].tolist():
            itrian.extend(trian[np.nonzero(trian == itri)[0]])
        itrian = np.array(itrian).astype('int')
        bool_itrian = np.in1d(itrian, ipos[0]).reshape(itrian.shape[0], 3)
        itrian[np.nonzero(bool_itrian == False)] = 0
        citri = np.vstack([np.vstack([itrian[:, 0], itrian[:, 1]]).T, np.vstack([itrian[:, 1], itrian[:, 2]]).T,
                           np.vstack([itrian[:, 2], itrian[:, 0]]).T])
        bcitri = (citri != 0).sum(1)
        valp = citri[bcitri == 2]
        # TODO: plot the picture and save them somewhere
        # fig = figure(figsize=(15, 15))
        # fig.suptitle('region ' + str(int(ilab)))
        # ax = fig.add_subplot(111, projection='3d')
        # xlabel('x')
        # ylabel('y')
        # for iv in np.arange(valp.shape[0]):
        #     ax.plot(vert[valp[iv], 0], vert[valp[iv], 1], vert[valp[iv], 2])
        # show()

    def _run_interface(self, runtime):

        vert = np.loadtxt(self.inputs.vertices)
        trian = np.loadtxt(self.inputs.triangles)
        texture = np.loadtxt(self.inputs.region_mapping)

        res = np.array(self.calculate_connected(texture, vert, trian))
        wrong_labels = res[np.where(res[:, 1] > 0.), 0][0]
        if len(wrong_labels) == 0:
            print 'everything is fine, continuing'
        else:
            print "WARNING, some region have several components"
            for iwrong in wrong_labels:
                # TODO: handle more than two components
                (V1, V2) = self.find_both_components(texture, vert, trian, iwrong)
                if self.inputs.check and len(self.inputs.display) > 0:
                    print "checking"
                    self.check_region_mapping(texture, vert, trian, iwrong)
                    i = 0
                    while True and i < 10:
                        try:
                            choice_user = raw_input("""Do you want to get rid of region with:
                                                1) {0} nodes
                                                2) {1} nodes
                                                3) continue the pipeline anyway
                                                (answer: 1, 2, or 3)? \n""".format(len(V1), len(V2)))
                            print "you chose " + choice_user
                            choice_user = int(choice_user)
                            if choice_user not in [1, 2, 3]:
                                raise ValueError
                            break
                        except ValueError:
                            print 'please choose 1, 2 or 3'
                            i += 1
                            continue
                    else:
                        print 'failure total, no check mode'
                else:
                    print "no check, selecting automatically the smallest components"
                    choice_user = np.argmin((len(V1), len(V2))) + 1
                if choice_user == 3:
                    print "keep that correction"
                    np.savetxt(self._get_outfilename(), texture)
                elif choice_user == 1:
                    texture = self.correct_sub_region(texture, trian, V1)
                elif choice_user == 2:
                    texture = self.correct_sub_region(texture, trian, V2)
                else:
                    print('failure of the choice')

            np.savetxt(self._get_outfilename(), texture)
        return runtime

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['region_mapping'] = self.inputs.region_mapping
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
            out_file_triangles = fname_presuffix(self.inputs.triangles[0], newpath=os.getcwd(),
                                                 suffix='_triangles', use_ext=True)
        if not isdefined(out_file_vertices):
            out_file_vertices = fname_presuffix(self.inputs.vertices[0], newpath=os.getcwd(),
                                                suffix='_vertices', use_ext=True)
        if not isdefined(out_file_texture):
            out_file_texture = fname_presuffix(self.inputs.textures[0], newpath=os.getcwd(),
                                               suffix='_texture', use_ext=True)
        return (os.path.abspath(out_file_triangles),
                os.path.abspath(out_file_vertices),
                os.path.abspath(out_file_texture))

    def _run_interface(self, runtime):
        (out_file_triangles, out_file_vertices, out_file_texture) = self._get_outfilename()
        lh_reg_map = np.loadtxt(self.inputs.textures[0])
        rh_reg_map = np.loadtxt(self.inputs.textures[1])
        lh_vert = np.loadtxt(self.inputs.vertices[0])
        rh_vert = np.loadtxt(self.inputs.vertices[1])
        lh_trian = np.loadtxt(self.inputs.triangles[0])
        rh_trian = np.loadtxt(self.inputs.triangles[1])
        vertices = np.vstack([lh_vert, rh_vert])
        triangles = np.vstack([lh_trian, rh_trian + lh_vert.shape[0]])
        region_mapping = np.hstack([lh_reg_map, rh_reg_map])
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
    label = traits.Str(desc="FreeSurfer Subject Id", argstr='-l %s',
                       exists=True, mandatory=True)
    subjects_dir = Directory(desc='path to subject directory', argstr='-d %s',
                             exists=True)
    freesurfer_home = Directory(desc='freesurfer home directory', argstr='-f %s',
                                exists=True)
    scripts_directory = Directory(desc='directory for scripts', exist=True)
    # out_file_subcortical_surf_list = traits.List(Files(desc="subcortical files"))


class Aseg2SrfOutputSpec(TraitedSpec):
    subcortical_surf_file = File(desc="Output subcortical surface", exists=True)


class Aseg2Srf(CommandLine):
    input_spec = Aseg2SrfInputSpec
    output_spec = Aseg2SrfOutputSpec
    label = "4"
    # "4, 5, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 26, 28, 43,\
    #  44, 46, 47, 49, 50, 51, 52, 53, 54, 58, 60, 251, 252, 253,\
    # 254, 255"
    _cmd = os.getcwd()

    def __init__(self, *args, **kwargs):
        super(Aseg2Srf, self).__init__(*args, **kwargs)
        self._cmd = self.get_scripts_dir() + "/aseg2srf"

    def get_scripts_dir(self):
        if isdefined(self.inputs.scripts_directory):
            return self.inputs.scripts_directory
        else:
            return os.getcwd()

    def _gen_subjects_dir(self):
        return os.getcwd()

    def _list_outputs(self):
        outputs = self.output_spec().get()
        # if isdefined(self.inputs.out_file_subcortical_surf_list):
        #    out_files = self.inputs.out_file_subcortical_surf_list
        # else:
        if isdefined(self.inputs.subjects_dir):
            sd = self.inputs.subjects_dir
        else:
            sd = os.environ('SUBJECTS_DIR')
        if len(self.inputs.label)==1:
            out_files = os.path.join(sd, self.inputs.subject_id, 'ascii', 'aseg_00' + self.inputs.label + '.srf')
        else:
            out_files = os.path.join(sd, self.inputs.subject_id, 'ascii', 'aseg_0' + self.inputs.label + '.srf')
        outputs['subcortical_surf_file'] = out_files
        return outputs


class ListSubcorticalInputSpect(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True)
    out_file_vertices_sub = traits.Str(desc="output subcortical surfaces")
    out_file_triangles_sub = traits.Str(desc="output subcortical triangles")


class ListSubcorticalOutputSpect(TraitedSpec):
    triangles_sub = File(exists=True)
    vertices_sub = File(exists=True)


class ListSubcortical(BaseInterface):
    input_spec = ListSubcorticalInputSpect
    output_spec = ListSubcorticalOutputSpect

    def _get_outfilename(self, in_file):
        out_file_vertices_sub = self.inputs.out_file_vertices_sub
        out_file_triangles_sub = self.inputs.out_file_triangles_sub
        if not isdefined(out_file_triangles_sub):
            out_file_triangles_sub = os.path.join(os.getcwd(), in_file +
                                                  '_triangles_sub.txt')
        if not isdefined(out_file_vertices_sub):
            out_file_vertices_sub = os.path.join(os.getcwd(), in_file +
                                                 '_vertices_sub.txt')
        return (os.path.abspath(out_file_triangles_sub),
                os.path.abspath(out_file_vertices_sub))

    def _run_interface(self, runtime):
        in_file = self.inputs.in_file
        f = open(in_file, 'r')
        f.readline()
        data = f.readline()
        g = data.split(' ')
        nb_vert = int(g[0])
        # nb_tri = int(g[1].split('\n')[0])
        f.close()
        a = np.loadtxt(in_file, skiprows=2, usecols=(0, 1, 2))
        vert = a[:nb_vert]
        tri = a[nb_vert:].astype('int')
        (out_file_triangles_sub, out_file_vertices_sub) = self._get_outfilename(in_file)
        np.savetxt(out_file_vertices_sub, vert)
        np.savetxt(out_file_triangles_sub, tri)
        return runtime

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['triangles_sub'] = self._get_outfilename(self.inputs.in_file)[0]
        outputs['vertices_sub'] = self._get_outfilename(self.inputs.in_file)[1]
        return outputs


class GrabDataCorInputSpec(BaseInterfaceInputSpec):
    subject_id = traits.Str(exists=True, mandatory=True)
    freesurfer_directory = Directory(exists=True, mandatory=True)
    out_file_aparcaseg = File(desc='output filename for aparcaseg file')
    out_file_t1 = File(desc='output filename for T1 file')


class GrabDataCorOutputSpec(TraitedSpec):
    aparcaseg = File(desc='output filename for aparcaseg file', exists=True)
    t1 = File(desc='output filename for t1', exists=True)


class GrabDataCor(BaseInterface):
    """
    grab the files produced by freesurfer.
    """
    input_spec = GrabDataCorInputSpec
    output_spec = GrabDataCorOutputSpec

    def _run_interface(self, runtime):
        return runtime

    def _list_outputs(self):
        path_fs_subj = os.path.join(os.path.abspath(self.inputs.freesurfer_directory), self.inputs.subject_id)
        out_file_aparcaseg = os.path.join(path_fs_subj, 'mri', 'aparc+aseg.mgz')
        out_file_t1 = os.path.join(path_fs_subj, 'mri', 'T1.mgz')
        outputs = self._outputs().get()
        outputs['aparcaseg'] = out_file_aparcaseg
        outputs['t1'] = out_file_t1
        return outputs


class ComputeConnectivityFilesInputSpec(BaseInterfaceInputSpec):
    verts = File(exists=True, mandatory=True, desc='surface vertices')
    tri = File(exists=True, mandatory=True, desc='surface triangles')
    region_mapping = File(exists=True, mandatory=True, desc='region mapping')
    weights = File(exists=True, mandatory=True, desc='weight matrix')
    tract_lengths = File(exists=True, mandatory=True, desc='tract length matrix')
    vertices_sub_list = traits.List(File(exists=True, mandatory=True))
    triangles_sub_list = traits.List(File(exists=True, mandatory=True))
    corr_table = traits.File(exists=True, mandatory=True, desc='correspondance table')
    name_regions = traits.File(exists=True, mandatory=True, desc='name regions')
    out_file_weights = File(desc='resulting weights')
    out_file_tract_lengths = File(desc='resulting tract lengths')
    out_file_areas = File(desc='resulting areas')
    out_file_average_orientations = File(desc='resulting average orientations')
    out_file_centres = File(desc='resulting centres')


class ComputeConnectivityFilesOutputSpec(TraitedSpec):
    weights = File(exists=True, desc='resulting weights')
    tract_lengths = File(exists=True, desc='resulting tract lengths')
    areas = File(exists=True, desc='resulting areas')
    average_orientations = File(exists=True, desc='resulting average orientations')
    centres = File(exists=True, desc='resulting centres')


class ComputeConnectivityFiles(BaseInterface):
    input_spec = ComputeConnectivityFilesInputSpec
    output_spec = ComputeConnectivityFilesOutputSpec

    def compute_triangle_areas(self, vertices, triangles):
        """Calculates the area of triangles making up a surface."""
        tri_u = vertices[triangles[:, 1], :] - vertices[triangles[:, 0], :]
        tri_v = vertices[triangles[:, 2], :] - vertices[triangles[:, 0], :]
        tri_norm = np.cross(tri_u, tri_v)
        triangle_areas = np.sqrt(np.sum(tri_norm ** 2, axis=1)) / 2.0
        triangle_areas = triangle_areas[:, np.newaxis]
        return triangle_areas

    def compute_region_areas(self, triangles_areas, vertex_triangles):
        avt = np.array(vertex_triangles)
        # NOTE: Slightly overestimates as it counts overlapping border triangles,
        #      but, not really a problem provided triangle-size << region-size.
        regs = map(set, avt)
        region_triangles = set.union(*regs)
        region_surface_area = triangles_areas[list(region_triangles)].sum()
        return region_surface_area

    def compute_region_orientation(self, vertex_normals):
        average_orientation = np.zeros((1, 3))
        # Average orientation of the region
        orient = vertex_normals[:, :]
        avg_orient = np.mean(orient, axis=0)
        average_orientation = avg_orient / np.sqrt(np.sum(avg_orient ** 2))
        region_orientation = average_orientation
        return region_orientation

    def compute_vertex_triangles(self, number_of_vertices, number_of_triangles, triangles):
        vertex_triangles = [[] for _ in xrange(number_of_vertices)]
        for k in xrange(number_of_triangles):
            vertex_triangles[triangles[k, 0]].append(k)
            vertex_triangles[triangles[k, 1]].append(k)
            vertex_triangles[triangles[k, 2]].append(k)
        return vertex_triangles

    def compute_vertex_normals(self, number_of_vertices, vertex_triangles, triangles,
                               triangle_angles, triangle_normals, vertices):
        """
        Estimates vertex normals, based on triangle normals weighted by the
        angle they subtend at each vertex...
        """
        vert_norms = np.zeros((number_of_vertices, 3))
        bad_normal_count = 0
        for k in xrange(number_of_vertices):
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

    def compute_triangle_angles(self, vertices, number_of_triangles, triangles):
        """
        Calculates the inner angles of all the triangles which make up a surface
        """
        verts = vertices
        # TODO: Should be possible with arrays, ie not nested loops...
        # A short profile indicates this function takes 95% of the time to compute normals
        # (this was a direct translation of some old matlab code)
        angles = np.zeros((number_of_triangles, 3))
        for tt in xrange(number_of_triangles):
            triangle = triangles[tt, :]
            for ta in xrange(3):
                ang = np.roll(triangle, -ta)
                angles[tt, ta] = np.arccos(np.dot(
                    (verts[ang[1], :] - verts[ang[0], :]) /
                    np.sqrt(np.sum((verts[ang[1], :] - verts[ang[0], :]) ** 2, axis=0)),
                    (verts[ang[2], :] - verts[ang[0], :]) /
                    np.sqrt(np.sum((verts[ang[2], :] - verts[ang[0], :]) ** 2, axis=0))))

        return angles

    def compute_triangle_normals(self, triangles, vertices):
        """Calculates triangle normals."""
        tri_u = vertices[triangles[:, 1], :] - vertices[triangles[:, 0], :]
        tri_v = vertices[triangles[:, 2], :] - vertices[triangles[:, 0], :]
        tri_norm = np.cross(tri_u, tri_v)

        try:
            triangle_normals = tri_norm / np.sqrt(np.sum(tri_norm ** 2, axis=1))[:, np.newaxis]
        except FloatingPointError:
            # TODO: NaN generation would stop execution, however for normals this case could maybe be
            # handled in a better way.
            triangle_normals = tri_norm
        return triangle_normals

    def compute_region_areas_cortex(self, triangle_areas, vertex_triangles, region_mapping):
        regions = np.unique(region_mapping)
        region_surface_area = np.zeros((np.max(np.unique(regions)) + 1, 1))
        avt = np.array(vertex_triangles)
        # NOTE: Slightly overestimates as it counts overlapping border triangles,
        #      but, not really a problem provided triangle-size << region-size.
        for k in regions:
            regs = map(set, avt[region_mapping == k])
            region_triangles = set.union(*regs)
            region_surface_area[k] = triangle_areas[list(region_triangles)].sum()

        return region_surface_area

    def compute_region_orientation_cortex(self, vertex_normals, region_mapping):
        regions = np.unique(region_mapping)
        average_orientation = np.zeros((np.max(np.unique(regions)) + 1, 3))
        # Average orientation of the region
        for k in regions:
            orient = vertex_normals[region_mapping == k, :]
            avg_orient = np.mean(orient, axis=0)
            average_orientation[k, :] = avg_orient / np.sqrt(np.sum(avg_orient ** 2))

        return average_orientation

    def compute_region_center_cortex(self, vertices, region_mapping):
        regions = np.unique(region_mapping)
        region_center = np.zeros((np.max(np.unique(region_mapping)) + 1, 3))
        # Average orientation of the region
        for k in regions:
            vert = vertices[region_mapping == k, :]
            region_center[k, :] = np.mean(vert, axis=0)

        return region_center

    def _get_outfilename(self):
        out_file_weights = self.inputs.out_file_weights
        out_file_tract_lengths = self.inputs.out_file_tract_lengths
        out_file_areas = self.inputs.out_file_areas
        out_file_average_orientations = self.inputs.out_file_average_orientations
        out_file_centres = self.inputs.out_file_centres

        if not isdefined(out_file_weights):
            out_file_weights = os.path.join(os.getcwd(), 'weights.txt')
        if not isdefined(out_file_tract_lengths):
            out_file_tract_lengths = os.path.join(os.getcwd(), 'tract_lengths.txt')
        if not isdefined(out_file_areas):
            out_file_areas = os.path.join(os.getcwd(), 'areas.txt')
        if not isdefined(out_file_average_orientations):
            out_file_average_orientations = os.path.join(os.getcwd(), 'average_orientations.txt')
        if not isdefined(out_file_centres):
            out_file_centres = os.path.join(os.getcwd(), 'centres.txt')

        return (os.path.abspath(out_file_weights),
                os.path.abspath(out_file_tract_lengths),
                os.path.abspath(out_file_areas),
                os.path.abspath(out_file_average_orientations),
                os.path.abspath(out_file_centres))

    def _run_interface(self, runtime):
        # Cortex
        (out_file_weights, out_file_tract_lengths, out_file_areas, out_file_average_orientations,
         out_file_centres) = self._get_outfilename()
        # import data
        verts = np.loadtxt(self.inputs.verts)
        tri = np.loadtxt(self.inputs.tri).astype(int)
        region_mapping = np.loadtxt(self.inputs.region_mapping)
        # save connectivity and tract length matrices
        weights = np.loadtxt(self.inputs.weights)
        tract_lengths = np.loadtxt(self.inputs.tract_lengths)
        weights = weights + weights.transpose() - np.diag(np.diag(weights))
        # add the first region
        weights = np.vstack([np.zeros((1, weights.shape[0])), weights])
        weights = np.hstack([np.zeros((weights.shape[0], 1)), weights])
        tract_lengths = tract_lengths + tract_lengths.transpose()  # because diagonal nul
        # tck2connectome produces the number of steps in the tracking file
        # so we have to multiply by the step size in mm which is 0.5
        tract_lengths *= 0.5
        tract_lengths = np.vstack([np.zeros((1, tract_lengths.shape[0])), tract_lengths])
        tract_lengths = np.hstack([np.zeros((tract_lengths.shape[0], 1)), tract_lengths])
        np.savetxt(out_file_weights, weights, fmt='%d')
        np.savetxt(out_file_tract_lengths, tract_lengths, fmt='%.3f')

        # compute centers
        centers = self.compute_region_center_cortex(verts, region_mapping)

        # calculate average orientations
        number_of_vertices = int(verts.shape[0])
        number_of_triangles = int(tri.shape[0])
        vertex_triangles = self.compute_vertex_triangles(number_of_vertices, number_of_triangles,
                                                         tri)
        triangle_normals = self.compute_triangle_normals(tri, verts)
        triangle_angles = self.compute_triangle_angles(verts, number_of_triangles, tri)
        vertex_normals = self.compute_vertex_normals(number_of_vertices, vertex_triangles,
                                                     tri, triangle_angles,
                                                     triangle_normals, verts)
        orientations = self.compute_region_orientation_cortex(vertex_normals, region_mapping)

        # compute areas
        triangle_areas = self.compute_triangle_areas(verts, tri)
        areas = self.compute_region_areas_cortex(triangle_areas, vertex_triangles, region_mapping)

        # subcorticals
        corr_table = np.loadtxt(self.inputs.corr_table)
        for sub_verts, sub_tri in zip(self.inputs.vertices_sub_list, self.inputs.triangles_sub_list):
            list_pos_val = ['16', '08', '10', '11', '12', '13', '17', '18', '26', '47', '49', '50', '51', '52', '53',
                            '54', '58']
            val = [iarg for iarg in list_pos_val if '_0' + iarg + '.' in sub_verts][0]
            verts = np.loadtxt(sub_verts)
            tri = np.loadtxt(sub_tri)
            tri = tri.astype(int)

            curr_center = np.mean(verts, axis=0)
            # import pdb; pdb.set_trace()
            indx = int(corr_table[np.nonzero(corr_table[:, 0] == np.int(val)), 1] - 1)
            centers[indx, :] = curr_center
            # Now calculate average orientations
            number_of_vertices = int(verts.shape[0])
            number_of_triangles = int(tri.shape[0])
            vertex_triangles = self.compute_vertex_triangles(number_of_vertices, number_of_triangles,
                                                             tri)
            triangle_normals = self.compute_triangle_normals(tri, verts)
            triangle_angles = self.compute_triangle_angles(verts, number_of_triangles, tri)
            vertex_normals = self.compute_vertex_normals(number_of_vertices, vertex_triangles,
                                                         tri, triangle_angles,
                                                         triangle_normals, verts)
            average_orientation = self.compute_region_orientation(vertex_normals)
            orientations[indx, :] = average_orientation

            triangle_areas = self.compute_triangle_areas(verts, tri)
            region_areas = self.compute_region_areas(triangle_areas, vertex_triangles)
            areas[indx] = region_areas

        # save orientations and areas
        np.savetxt(out_file_areas, areas, fmt='%.2f')
        np.savetxt(out_file_average_orientations,
                   orientations, fmt='%.2f %.2f %.2f')

        # add the name to centers
        f = open(self.inputs.name_regions, 'rb')
        list_name = []
        for line in f:
            list_name.append(line)
        f.close()

        f = open(out_file_centres, 'w')
        for i, name in enumerate(list_name):
            f.write(str(name[:-1]) + ' ')
            for j in range(3):
                f.write('{:.4f} '.format(centers[i, j]))
            f.write('\n')
        f.close()

        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        (out_file_weights,
         out_file_tract_lengths,
         out_file_areas,
         out_file_average_orientations,
         out_file_centres) = self._get_outfilename()
        outputs['weights'] = out_file_weights
        outputs['tract_lengths'] = out_file_tract_lengths
        outputs['areas'] = out_file_areas
        outputs['average_orientations'] = out_file_average_orientations
        outputs['centres'] = out_file_centres
        return outputs
