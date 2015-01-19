import nipype.interfaces.fsl as fsl
import nipype.interfaces.mrtrix as mrt
import nipype.interfaces.freesurfer as fs
import nipype.pipeline.engine as pe
from nipype.workflows.dmri.fsl.artifacts import ecc_pipeline 
import scripts.utility as su 

def connecitivity_workflow():
    convert_dicom2nii = pe.Node(interface=mrt.MRConvert, name='convert_dicom2nii')
    extract_bvecs_bvals = pe.Node(interface=mtr3.utils.MRInfo, name='extract_bvecs_bvals')
# ecc = ecc_pipeline()
    convert_nii2dwi = pe.Node(interface=mrt.MRConvert, name='convert_nii2dwi')
    create_mask = pe.Node(interface=mrt3.dwi2mask, name='create_mask')
    dwiextract = .Node(interface=mrt3.dwiextract, name='dwiextract')
    dwi2nii = pe.Node(interface=mrt.MRConvert, name='dwi2nii')
    T1_mgz2nii = pe.Node(interface=fsl.MRIConvert, name='T1_mgz2nii')
    aparcaseg_mgz2nii = pe.Node(interface=fs.preprocess.MRIConvert, name='aparcaseg_mgz2nii')
    reorient2std = pe.Node(interface=fsl.utils.Reorient2Std, name='reorient2std')
    diff2struct = pe.Node(interface=fsl.preprocess.FLIRT, name='diff2struct')
    convertxfm = pe.Node(interface=fsl.utils.ConvertXF, name='convertXFM')
    aparcaseg2diff = pe.Node(interface=fsl.preprocess.ApplyXFM, name='aparcaseg2diff')
    T12diff = pe.Node(interface=fsl.preprocess.ApplyXFM, name='T12diff')
    T12diff.inputs.interp = 'nearestneighbour'
    dwi2response = pe.Node(interface=mrt3.preprocess.Dwi2Response, name='dwi2response')
    dwi2fod = pe.Node(interface=mrt3.preprocess.Dwi2Fod, name='dwi2fod')
    dwi2fod.inputs.lmax = 8
    act_anat_prepare_fls = pe.node(interface=mrt3.ActAnatPrepareFSL, name='act_anat_prepare_fsl')
    5tt2gmwmi = pe.Node(interface=mrt3.5tt2Gmwmi, name='5tt2gmwmi')
    tckgen = pe.node(interface=mrt3.tracking.Tckgen, name='tckgen')
    tckgen.inputs.unidirectional = True
    tckgen.inputs.seed_gmwmi = iFOD2
    tckgen.inputs.maxlength = 150
    tckgen.inputs.num = 5000000
    tcksift = pe.Node(interface=mrt3.TckSift, neame='tcksift')
    labelconfig = pe.Node(interface=mrt3.utils.LabelConfig, name='labelconfig')
    labelconfig.inputs.lut_freesurfer = 
    tck2connectome_weights = pe.Node(interface=mrt3.tracking.Tck2Connectome, name='tck2connectome_weights')
    tck2connectome_tract_lengths = pe.Node(interface=mrt3.tracking.Tck2Connectome, name='tck2connectome_tract_lengths')
    compute_connectivity = pe.Node(interface=su.ComputeConnectivityFiles, name='compute_connecitivty')







