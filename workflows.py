import nipype.interfaces.fsl as fsl
import nipype.interfaces.mrtrix as mrt
import nipype.pipeline.engine as pe
from nipype.workflows.dmri.fsl.artifacts import ecc_pipeline 

pe.Node(interface=mrt.MRConvert, name='convert_dicom2nii')
pe.Node(interface=mtr3.mrinfo, name='extract_bvec_bvals')
ecc = ecc_pipeline()
pe.Node(interface=mrt.MRConvert, name='convert_nii2dwi')
pe.Node(interface=mrt3.dwi2mask, name='create_mask')
pe.Node(interface=mrt3.dwiextract, name='dwiextract')
pr.Node(interface=mrt.mrconvert, name='dwi2nii')




