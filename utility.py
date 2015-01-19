class Asegg2SrfInputSpec(CommandLineSpec):
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

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['output_volume'] = os.path.abspath(self.inputs.in_subject_id  + '/scaii/seg_.asc')
        return outputs
