import nipype.interfaces.fsl as fsl
import nipype.interfaces.mrtrix as mrt
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as niu
import nipype.pipeline.engine as pe
from nipype.workflows.dmri.fsl.artifacts import ecc_pipeline 
import scripts.utility as su 

def DiffusionTractography(name="connectivity"):
    inputnode = pe.Node(niu.IdentityInterface(fields=['in_file']), name='inputnode')
    convert_dicom2nii = pe.Node(interface=mrt.MRConvert, name='convert_dicom2nii')
    extract_bvecs_bvals = pe.Node(interface=mtr3.utils.MRInfo, name='extract_bvecs_bvals')
# ecc = ecc_pipeline()
    convert_nii2dwi = pe.Node(interface=mrt.MRConvert, name='convert_nii2dwi')
    create_mask = pe.Node(interface=mrt3.Dwi2Mask, name='create_mask')
    dwi_extract_lowb = .Node(interface=mrt3.DwiExtract, name='dwi_extract_lowb')
    lowb_mif2lowb_nii = pe.Node(interface=mrt.MRConvert, name='lowb_mif2lowb_nii')
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

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode, convert_dicom2nii, [('in_file', 'in_file')]),
        (convert_dicom2nii, extract_bvecs_bvals, [('converted', 'in_files')]),
        (convert_dicom2nii, convert_nii2dwi, [('converted', 'in_file'),
                                              ('bvecs', 'bvals', 'fslgrad')]),
        (convert_nii2dwi, create_mask, [('converted', 'in_file')]),
        (convert_nii2dwi, dwi_extract_lowb, [('converted', 'in_file')]),
        (create_mask, dwi_extract_lowb, [('out_file', 'in_file')]),
        (dwi_extract_lowb, lowb_mif2lowb_nii, [('out_file', 'in_file')]),
        (

def Coregistration(name='coregistration'):
    inputnode = pe.Node(niu.IdentityInterface(fields=['in_T1.mgz', 'in_lowb.nii',
                        'in_aparcaseg.mgz']), name='inputnode')
    T1_mgz2nii = pe.Node(interface=fs.preprocess.MRIConvert, name='T1_mgz2nii')
    T1_mgz2nii.inputs.in_type = 'mgz'
    T1_mgz2nii.inputs.out_type = 'nii'
    T1_mgz2nii.inputs.out_orientation = 'RAS'
    aparcaseg_mgz2nii = pe.Node(interface=fs.preprocess.MRIConvert, name='aparcaseg_mgz2nii')
    aparcaseg_mgz2nii.inputs.in_type = 'mgz'
    aparcaseg_mgz2nii.inputs.out_type = 'nii'
    aparcaseg_mgz2nii.inputs.out_orientation = 'RAS'
    reorient2std = pe.Node(interface=fsl.utils.Reorient2Std, name='reorient2std')
    diff2struct = pe.Node(interface=fsl.preprocess.FLIRT, name='diff2struct')
    diff2struct.inputs.searchr_x = [180, 180]
    diff2struct.inputs.searchr_y = [180, 180]
    diff2struct.inputs.searchr_z = [180, 180]
    diff2struct.inputs.cost = 'mutualinfo'
    convertxfm = pe.Node(interface=fsl.utils.ConvertXFM, name='convertXFM')
    convertxfm.inputs.invert_xfm = True
    aparcaseg2diff = pe.Node(interface=fsl.preprocess.ApplyXFM, name='aparcaseg2diff')
    aparcaseg2diff.inputs.apply_xfm = True
    aparcaseg2diff.inputs.interp = 'nearestneighbour'
    T12diff = pe.Node(interface=fsl.preprocess.ApplyXFM, name='T12diff')
    T12diff.inputs.apply_xfm = True
    T12diff.inputs.interp = 'nearestneighbour'
    outputnode = pe.Node(niu.IdentityInterface(fields=['out_file'], name='outputnode'))
    
    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode, T1_mgz2nii, [('in_T1.mgz', 'in_file')]),
        (inputnode, aparcaseg_mgz2nii, [('in_aparcaseg.mgz', 'in_file')]),
        (aparcaseg_mgz2nii, reorient2std, [('out_file', 'in_file')]),
        (inputnode, diff2struct, [('in_lowb.nii', 'in_file')]),
        (T1_mgz2nii, diff2struct, [('out_file', 'reference')]),
        (diff2struct, convertxfm, [('out_matrix_file', 'in_file')]),
        (aparcaseg_mgz2nii, aparcaseg2diff, [('out_file', 'in_file')]),
        (inputnode, aparcaseg2diff, [('in_lowb.nii', 'reference')]),
        (convertxfm, aparcaseg2diff, [('out_file', 'in_matrix_file')]),
        (T1_mgz2nii, T12diff, [('out_file', 'in_file')]),
        (inputnode, T12diff, [('in_lowb.nii', 'reference')]),
        (convertxfm, T12diff, [('out_file', 'in_matrix_file')]),

    return wf







