from nipype.interfaces.base import (TraitedSpec, File, traits, CommandLineInputSpec, CommandLine, isdefined,
                                    InputMultiPath)
from nipype.utils.filemanip import split_filename
import os
import os.path as op


class MRConvertInputSpec(CommandLineInputSpec):
    in_file = File(exists=True, argstr='%s', mandatory=True, position=-2,
                   desc='voxel-order data filename')
    out_filename = File(genfile=True, argstr='%s', position=-1, desc='Output filename')
    extract_at_axis = traits.Enum(1, 2, 3, argstr='-coord %s', position=1,
                                  desc='Extract data only at the coordinates specified. This option specifies the'
                                       ' Axis. Must be used in conjunction with extract_at_coordinate.')
    extract_at_coordinate = traits.List(traits.Float, argstr='%s', sep=',', position=2, minlen=1, maxlen=3,
                                        desc='"Extract data only at the coordinates specified. This option specifies '
                                             'the coordinates. Must be used in conjunction with extract_at_axis. Three '
                                             'comma-separated numbers giving the size of each voxel in mm.')
    voxel_dims = traits.List(traits.Float, argstr='-vox %s', sep=',',
                             position=3, minlen=3, maxlen=3,
                             desc='Three comma-separated numbers giving the size of each voxel in mm.')
    output_datatype = traits.Enum("nii", "float", "char", "short", "int", "long", "double", argstr='-output %s',
                                  position=2,
                                  desc='"i.e. Bfloat". Can be "char", "short", "int", "long", "float" or "double"')
    # , usedefault=True)
    extension = traits.Enum("mif", "nii", "float", "char", "short", "int", "long", "double", position=2,
                            desc='"i.e. Bfloat". Can be "char", "short", "int", "long", "float" or "double"',
                            usedefault=True)
    layout = traits.Enum("nii", "float", "char", "short", "int", "long", "double", argstr='-output %s', position=2,
                         desc='specify the layout of the data in memory. The actual layout produced will depend on '
                              'whether the output image format can support it.')
    resample = traits.Float(argstr='-scale %d', position=3,
                            units='mm', desc='Apply scaling to the intensity values.')
    offset_bias = traits.Float(argstr='-scale %d', position=3,
                               units='mm', desc='Apply offset to the intensity values.')
    replace_NaN_with_zero = traits.Bool(argstr='-zero', position=3, desc="Replace all NaN values with zero.")
    prs = traits.Bool(argstr='-prs', position=3,
                      desc="Assume that the DW gradients are specified in the PRS frame (Siemens DICOM only).")
    grad = traits.Str(argstr='-grad %s', position=3, desc='mrtrix encoding')


class MRConvertOutputSpec(TraitedSpec):
    converted = File(exists=True, desc='path/name of 4D volume in voxel order')


class MRConvert(CommandLine):
    """
    Perform conversion between different file types and optionally extract a subset of the input image.

    If used correctly, this program can be a very useful workhorse.
    In addition to converting images between different formats, it can
    be used to extract specific studies from a data set, extract a specific
    region of interest, flip the images, or to scale the intensity of the images.

    Example
    -------

    >>> import nipype.interfaces.mrtrix as mrt
    >>> mrconvert = mrt.MRConvert()
    >>> mrconvert.inputs.in_file = 'dwi_FA.mif'
    >>> mrconvert.inputs.out_filename = 'dwi_FA.nii'
    >>> mrconvert.run()                                 # doctest: +SKIP
    """

    _cmd = 'mrconvert'
    input_spec = MRConvertInputSpec
    output_spec = MRConvertOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['converted'] = self.inputs.out_filename
        if not isdefined(outputs['converted']):
            outputs['converted'] = op.abspath(self._gen_outfilename())
        else:
            outputs['converted'] = op.abspath(outputs['converted'])
        return outputs

    def _gen_filename(self, name):
        if name is 'out_filename':
            return self._gen_outfilename()
        else:
            return None

    def _gen_outfilename(self):
        _, name, _ = split_filename(self.inputs.in_file)
        if isdefined(self.inputs.out_filename):
            outname = self.inputs.out_filename
        else:
            outname = name + '_mrconvert.' + self.inputs.extension
        return outname


class DWI2TensorInputSpec(CommandLineInputSpec):
    in_file = InputMultiPath(File(exists=True), argstr='%s', mandatory=True,
                             position=-2, desc='Diffusion-weighted images')
    out_filename = File(name_template="%s_tensor.mif", name_source="in_file",
                        output_name="tensor", argstr='%s',
                        desc='Output tensor filename', position=-1)
    encoding_file = File(argstr='-grad %s', position=2,
                         desc=('Encoding file supplied as a 4xN text file with '
                               'each line is in the format [ X Y Z b ], where '
                               '[ X Y Z ] describe the direction of the applied '
                               'gradient, and b gives the b-value in units '
                               '(1000 s/mm^2). See FSL2MRTrix()'))
    ignore_slice_by_volume = traits.List(traits.Int, argstr='-ignoreslices %s',
                                         sep=' ', position=2, minlen=2,
                                         maxlen=2,
                                         desc=('Requires two values (i.e. [34 '
                                               '1] for [Slice Volume] Ignores '
                                               'the image slices specified '
                                               'when computing the tensor. '
                                               'Slice here means the z '
                                               'coordinate of the slice to be '
                                               'ignored.'))
    ignore_volumes = traits.List(traits.Int, argstr='-ignorevolumes %s',
                                 sep=' ', position=2, minlen=1,
                                 desc=('Requires two values (i.e. [2 5 6] for '
                                       '[Volumes] Ignores the image volumes '
                                       'specified when computing the tensor.'))
    quiet = traits.Bool(argstr='-quiet', position=1,
                        desc=("Do not display information messages or progress "
                              "status."))
    debug = traits.Bool(argstr='-debug', position=1,
                        desc="Display debugging messages.")


class DwiExtractInputSpec(CommandLineInputSpec):
    in_file = File(exists=True, mandatory=True, position=-2, argstr='%s', desc='the input DW image.')
    out_file = File(exists=True, genfile=True, position=-1, argstr='%s', desc='the output image (diffusion-weighted'
                                                                              ' volumes by default.')
    bzero = traits.Bool(argstr='-bzero', desc='output b=0 volumes instead of the diffusion weighted volumes.')


class DwiExtractOutputSpec(TraitedSpec):
    out_file = File(desc='the output image (diffusion-weighted'
                         ' volumes by default.', exists=True)


class DwiExtract(CommandLine):
    input_spec = DwiExtractInputSpec
    output_spec = DwiExtractOutputSpec
    _cmd = 'dwiextract'

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
            return name + '_dwiextracted' + ext

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_file'] = os.path.abspath(self._gen_outfilename())
        return outputs


class Fivett2GmwmiInputSpec(CommandLineInputSpec):
    fivett_in = File(exists=True, mandatory=True, position=-2, argstr='%s',
                     desc='the input 5TT segmented anatomical image')
    mask_out = File(exists=True, genfile=True, position=-1, argstr='%s', desc='the output mask image')


class Fivett2GmwmiOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='the output mask image')


class Fivett2Gmwmi(CommandLine):
    input_spec = Fivett2GmwmiInputSpec
    output_spec = Fivett2GmwmiOutputSpec
    _cmd = '5tt2gmwmi'

    def _gen_filename(self, name):
        if name is 'mask_out':
            return self._gen_outfilename()
        else:
            return None

    def _gen_outfilename(self):
        if isdefined(self.inputs.mask_out):
            return self.inputs.mask_out
        else:
            _, name, ext = split_filename(self.inputs.fivett_in)
            return name + '_gmwmi' + ext

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_file'] = os.path.abspath(self._gen_outfilename())
        return outputs


class TckSiftInputSpec(CommandLineInputSpec):
    in_tracks = File(exists=True, mandatory=True, position=-3, argstr='%s', desc='the input track file.')
    in_fod = File(exists=True, mandatory=True, position=-2, argstr='%s',
                  desc='input image containing the spherical harmonics of the fibre orientation distributions')
    out_tracks = File(exists=True, genfile=True, position=-1, argstr='%s', desc='the output filtered tracks file.')
    term_number = traits.Int(argstr='-term_number %s',
                             desc='number of streamlines - continue filtering until this number of'
                                  'streamlines remain')


class TckSiftOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='the output filtered tracks file.')


class TckSift(CommandLine):
    input_spec = TckSiftInputSpec
    output_spec = TckSiftOutputSpec
    _cmd = 'tcksift'

    def _gen_filename(self, name):
        if name is 'out_tracks':
            return self._gen_outfilename()
        else:
            return None

    def _gen_outfilename(self):
        if isdefined(self.inputs.out_tracks):
            return self.inputs.out_tracks
        else:
            _, name, ext = split_filename(self.inputs.in_tracks)
            return name + '_tcksifted' + ext

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_file'] = os.path.abspath(self._gen_outfilename())
        return outputs
