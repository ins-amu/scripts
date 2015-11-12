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


class DwiExtractInputSpec(CommandLineInputSpec):
    in_file = File(exists=True, mandatory=True, position=-2, argstr='%s', desc='the input DW image.')
    out_file = File(exists=True, genfile=True, position=-1, argstr='%s', desc='the output image (diffusion-weighted'
                                                                              ' volumes by default.')
    bzero = traits.Bool(argstr='bzero', desc='output b=0 volumes instead of the diffusion weighted volumes.')


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


class Fivett2GmwmiInputSpec():
    fivett_in = File(exists=True, mandatory=True, position=-2, argstr='%s',
                     desc='the input 5TT segmented anatomical image')
    mask_out = File(exists=True, genfile=True, position=-1, argstr='%s', desc='the output mask image')


class Fivett2GmwmiOutputSpec():
    out_file = File(exists=True, desc='the output mask image')


class Fivett2Gmwmi():
    input_spec = Fivett2GmwmiInputSpec
    output_spec = Fivett2GmwmiOutputSpec
    _cmd = '5tt2gmwmi'

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
            return name + '_tcksifted' + ext

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_file'] = os.path.abspath(self._gen_outfilename())
        return outputs
