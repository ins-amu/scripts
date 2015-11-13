import nipype.interfaces.fsl as fsl
import nipype.interfaces.mrtrix as mrt
import nipype.interfaces.mrtrix3 as mrt3
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as niu
import nipype.pipeline.engine as pe
from nipype.workflows.dmri.fsl.artifacts import ecc_pipeline
import utility as su
import mrtrix3_utility as mrt3u


def reconall():
    pass


def surface(name='surface'):
    """
    Surface workflow

    Generate vertices.txt, triangles.txt and region_mapping.txt
    """

    inputnode = pe.Node(interface=niu.IdentityInterface(fields=['pials', 'annots', 'ref_tables']), name='inputnode')
    # inputnode.iterables = [('pials', ["/Users/timp/Work/data/af/surf/rh.pial",
    #                                  "/Users/timp/Work/data/af//surf/lh.pial"]),
    #                       ('annots', ["/Users/timp/Work/data/af/label/rh.aparc.annot",
    #                                   "/Users/timp/Work/data/af/label/lh.aparc.annot"]),
    #                       ('ref_tables', ["/Users/timp/Desktop/scripts/rh_ref_table.txt",
    #                                       "/Users/timp/Desktop/scripts/lh_ref_table.txt"])]
    inputnode.synchronize = True

    pial2asc = pe.Node(interface=su.MRIsConvert(), name='pial2asc')
    pial2asc.inputs.out_datatype = 'asc'
    extract_high = pe.Node(interface=su.Fs2Txt(), name='extract_high')
    txt2off = pe.Node(interface=su.Txt2Off(), name='txt2off')
    remesher = pe.Node(interface=su.Remesher(), name='remesher')
    off2txt = pe.Node(interface=su.Off2Txt(), name='off2txt')
    region_mapping = pe.Node(interface=su.RegionMapping(), name='region_mapping')
    # region_mapping.inputs.scripts_directory = '/Users/timp/Desktop/scripts/'
    correct_region_mapping = pe.Node(interface=su.CorrectRegionMapping(), name='correct_region_mapping')
    check_region_mapping = pe.Node(interface=su.CheckRegionMapping(), name='check_region_mapping')
    # check_region_mapping.inputs.check = False
    # check_region_mapping.inputs.display = False
    # check_region_mapping.inputs.scripts_directory = '/Users/timp/Desktop/scripts/'
    reunify_both_hemi = pe.JoinNode(interface=su.ReunifyBothHemisphere(), joinsource="inputnode",
                                    joinfield=['vertices', 'triangles', 'textures'], name='reunify_both_hemi')
    output_node = pe.Node(interface=niu.IdentityInterface(fields=['texture', 'vertices', 'triangles']),
                          name='outputnode')

    # Connect workflow
    wf = pe.Workflow(name='workflow_surface')
    wf.connect([
        (inputnode, pial2asc, [('pials', 'in_file')]),
        (inputnode, region_mapping, [('annots', 'aparc_annot'),
                                     ('ref_tables', 'ref_table')]),
        (pial2asc, extract_high, [('converted', 'surface')]),
        (extract_high, txt2off, [('vertices', 'vertices'),
                                 ('triangles', 'triangles')]),
        (extract_high, region_mapping, [('vertices', 'vertices')]),
        (txt2off, remesher, [('surface_off', 'in_file')]),
        (remesher, off2txt, [('out_file', 'surface_off')]),
        # work to do on region mapping
        (off2txt, region_mapping, [('vertices_txt', 'vertices_downsampled'),
                                   ('triangles_txt', 'triangles_downsampled')]),
        (off2txt, correct_region_mapping, [('vertices_txt', 'vertices'),
                                           ('triangles_txt', 'triangles')]),
        (off2txt, check_region_mapping, [('vertices_txt', 'vertices'),
                                         ('triangles_txt', 'triangles')]),
        (off2txt, reunify_both_hemi, [('vertices_txt', 'vertices'),
                                      ('triangles_txt', 'triangles')]),
        (region_mapping, correct_region_mapping, [('out_file', 'texture')]),
        (correct_region_mapping, check_region_mapping, [('texture_corrected', 'region_mapping')]),
        (check_region_mapping, reunify_both_hemi, [('region_mapping', 'textures')]),
        (reunify_both_hemi, output_node, [('texture', 'region_mapping'),
                                          ('vertices', 'vertices'),
                                          ('triangles', 'triangles')])
    ])

    return wf


def subcorticalsurface(name="subcorticalsurfaces"):
    """
    Extraction of the subcortical surfaces from FreeSurfer

    Output a list of subcortical surfaces
    """
    inputnode = pe.Node(interface=niu.IdentityInterface(fields=['subject_id', 'labels', 'subjects_dir']),
                        name='inputnode')
    inputnode.inputs.subject_id = '151526'
    inputnode.inputs.subjects_dir = '/Users/timp/Desktop/scripts/tests/test_data/freesurfer_dir/'
    inputnode.iterables = [('labels', ['16', '08', '10', '11', '12', '13', '17', '18', '26', '47', '49', '50', '51',
                                       '52', '53', '54', '58'])]
    # inputnode.iterables = [('labels', ['16'])]
    aseg2srf = pe.Node(interface=su.Aseg2Srf(), name='aseg2srf')
    list_subcortical = pe.Node(interface=su.ListSubcortical(), name='list_subcortical')
    outputnode = pe.JoinNode(interface=niu.IdentityInterface(fields=['triangles_sub', 'vertices_sub']),
                             name='outputnode', joinsource="inputnode", joinfield=['triangles_sub', 'vertices_sub'])
    # Connect workflow
    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode, aseg2srf, [('subject_id', 'subject_id'),
                               ('labels', 'label'),
                               ('subjects_dir', 'subjects_dir')]),
        (aseg2srf, list_subcortical, [('subcortical_surf_file', 'in_file')]),
        (list_subcortical, outputnode, [('triangles_sub', 'triangles_sub'),
                                        ('vertices_sub', 'vertices_sub')]),
    ])
    return wf


# def Connectivity(name="connectivity"):
#     inputnode = pe.Node(interface=niu.IdentityInterface(fields=['in_file']), name='inputnode')
#     # TODO: first import
#     convert_dicom2nii = pe.Node(interface=mrt.MRConvert(), name='convert_dicom2nii')
#     extract_bvecs_bvals = pe.Node(interface=mtr3.utils.MRInfo(), name='extract_bvecs_bvals')
#     # ecc = ecc_pipeline()
#     convert_nii2dwi = pe.Node(interface=mrt.MRConvert(), name='convert_nii2dwi')
#     # END TODO

def tractography(name="tractography"):
    # subcort = subcorticalsurface()
    inputnode = pe.Node(interface=niu.IdentityInterface(
        fields=['converted', 'verts', 'tri', 'region_mapping', 'vertices_sub_list', 'triangles_sub_list', 'in_T1',
                'in_aparcaseg', 'in_lowb']),
        name='inputnode')
    # inputnode.inputs.converted = "/Users/timp/Work/data/af/connectivity/dwi.mif"
    # inputnode.inputs.verts = "/Users/timp/Work/data/af/af/surface/vertices.txt"
    # inputnode.inputs.tri = "/Users/timp/Work/data/af/af/surface/triangles.txt"
    # inputnode.inputs.region_mapping = "/Users/timp/Work/data/af/af/surface/region_mapping.txt"
    # inputnode.inputs.vertices_sub_list = "/Users/timp/Work/data/af/connectivity/dwi.mif"
    # inputnode.inputs.triangles_sub_list = "/Users/timp/Work/data/af/connectivity/dwi.mif"
    # inputnode.inputs.in_T1 = "/Users/timp/Work/data/af/mri/T1.mgz"
    # inputnode.inputs.in_aparcaseg = "/Users/timp/Work/data/af/mri/aparc+aseg.mgz"
    create_mask = pe.Node(interface=mrt3.utils.BrainMask(), name='create_mask')
    dwi_extract_lowb = pe.Node(interface=mrt3u.DwiExtract(), name='dwi_extract_lowb')
    dwi_extract_lowb.inputs.bzero = True
    lowb_mif2lowb_nii = pe.Node(interface=mrt.preprocess.MRConvert(), name='lowb_mif2lowb_nii')
    lowb_mif2lowb_nii.inputs.out_filename = 'lowb.nii'
    cor = Coregistration()
    dwi2response = pe.Node(interface=mrt3.preprocess.ResponseSD(), name='dwi2response')
    dwi2fod = pe.Node(interface=mrt3.reconst.EstimateFOD(), name='dwi2fod')
    dwi2fod.inputs.max_sh = 8
    act_anat_prepare_fsl = pe.Node(interface=mrt3.preprocess.ACTPrepareFSL(), name='act_anat_prepare_fsl')
    fivett2gmwmi = pe.Node(interface=mrt3u.Fivett2Gmwmi(), name='5tt2gmwmi')
    tckgen = pe.Node(interface=mrt3.tracking.Tractography(), name='tckgen')
    tckgen.inputs.unidirectional = True
    tckgen.inputs.algorithm = 'iFOD2'
    tckgen.inputs.max_length = 250.
    tckgen.inputs.step_size = 0.5
    # tckgen.inputs.n_tracks = 5000000
    tckgen.inputs.n_tracks = 5000
    tcksift = pe.Node(interface=mrt3u.TckSift(), name='tcksift')
    # tcksift.inputs.term_number = 2500000
    tcksift.inputs.term_number = 2500
    labelconfig = pe.Node(interface=mrt3.connectivity.LabelConfig(), name='labelconfig')
    labelconfig.inputs.in_config = '/Users/timp/Desktop/scripts/fs_region.txt'
    labelconfig.inputs.lut_fs = '/Applications/freesurfer/FreeSurferColorLUT.txt'
    tck2connectome_weights = pe.Node(interface=mrt3.connectivity.BuildConnectome(), name='tck2connectome_weights')
    tck2connectome_weights.inputs.search_radius = 2.
    tck2connectome_tract_lengths = pe.Node(interface=mrt3.connectivity.BuildConnectome(),
                                           name='tck2connectome_tract_lengths')
    tck2connectome_tract_lengths.inputs.search_radius = 2
    tck2connectome_tract_lengths.inputs.zero_diagonal = True
    tck2connectome_tract_lengths.inputs.metric = 'meanlength'
    compute_connectivity = pe.Node(interface=su.ComputeConnectivityFiles(), name='compute_connectivity')
    compute_connectivity.inputs.corr_table = '/Users/timp/Desktop/scripts/correspondance_mat.txt'
    compute_connectivity.inputs.name_regions = '/Users/timp/Desktop/scripts/name_regions.txt'
    output_node = pe.Node(
        interface=niu.IdentityInterface(
            fields=(['weights', 'tract_lengths', 'areas', 'centres', 'average_orientations'])), name='outputnode')

    wf = pe.Workflow(name=name)
    wf.connect([
        # (inputnode, convert_dicom2nii, [('in_file', 'in_file')]),
        # (convert_dicom2nii, extract_bvecs_bvals, [('converted', 'in_files')]),
        # (convert_dicom2nii, convert_nii2dwi, [('converted', 'in_file'),
        #                                      ('bvecs', 'bvals', 'fslgrad')]),
        # (convert_nii2dwi, create_mask, [('converted', 'in_file')]),
        # (subcort, inputnode, [('outputnode.triangles_sub', 'triangles_sub_list')]),
        # (subcort, inputnode, [('outputnode.vertices_sub', 'vertices_sub_list')]),
        (inputnode, create_mask, [('converted', 'in_file')]),
        (inputnode, dwi_extract_lowb, [('converted', 'in_file')]),
        (dwi_extract_lowb, lowb_mif2lowb_nii, [('out_file', 'in_file')]),
        (inputnode, cor, [('in_T1', 'inputnode.in_T1')]),
        (inputnode, cor, [('in_aparcaseg', 'inputnode.in_aparcaseg')]),
        (lowb_mif2lowb_nii, cor, [('converted', 'inputnode.in_lowb')]),
        (inputnode, dwi2response, [('converted', 'in_file')]),
        (create_mask, dwi2response, [('out_file', 'in_mask')]),
        (inputnode, dwi2fod, [('converted', 'in_file')]),
        (dwi2response, dwi2fod, [('out_file', 'response')]),
        (create_mask, dwi2fod, [('out_file', 'in_mask')]),
        (cor, act_anat_prepare_fsl, [('outputnode.out_T1_diff', 'in_file')]),
        (act_anat_prepare_fsl, fivett2gmwmi, [('out_file', 'fivett_in')]),
        (dwi2fod, tckgen, [('out_file', 'in_file')]),
        (fivett2gmwmi, tckgen, [('out_file', 'seed_gmwmi')]),
        (act_anat_prepare_fsl, tckgen, [('out_file', 'act_file')]),
        (tckgen, tcksift, [('out_file', 'in_tracks')]),
        (dwi2fod, tcksift, [('out_file', 'in_fod')]),
        (cor, labelconfig, [('outputnode.out_aparcaseg_diff', 'in_file')]),
        (tcksift, tck2connectome_weights, [('out_file', 'in_file')]),
        (labelconfig, tck2connectome_weights, [('out_file', 'in_parc')]),
        (tcksift, tck2connectome_tract_lengths, [('out_file', 'in_file')]),
        (labelconfig, tck2connectome_tract_lengths, [('out_file', 'in_parc')]),
        (tck2connectome_weights, compute_connectivity, [('out_file', 'weights')]),
        (tck2connectome_tract_lengths, compute_connectivity, [('out_file', 'tract_lengths')]),
        (inputnode, compute_connectivity, [('verts', 'verts')]),
        (inputnode, compute_connectivity, [('tri', 'tri')]),
        (inputnode, compute_connectivity, [('region_mapping', 'region_mapping')]),
        (inputnode, compute_connectivity, [('vertices_sub_list', 'vertices_sub_list')]),
        (inputnode, compute_connectivity, [('triangles_sub_list', 'triangles_sub_list')]),
        (compute_connectivity, output_node, [('weights', 'weights'),
                                             ('tract_lengths', 'tract_lengths'),
                                             ('areas', 'areas'),
                                             ('centres', 'centres'),
                                             ('average_orientations', 'average_orientations')])
    ])

    return wf


def Coregistration(name='coregistration'):
    inputnode = pe.Node(niu.IdentityInterface(fields=['in_T1', 'in_lowb',
                                                      'in_aparcaseg']), name='inputnode')
    T1_mgz2nii = pe.Node(interface=fs.preprocess.MRIConvert(), name='T1_mgz2nii')
    T1_mgz2nii.inputs.in_type = 'mgz'
    T1_mgz2nii.inputs.out_type = 'nii'
    T1_mgz2nii.inputs.out_orientation = 'RAS'
    aparcaseg_mgz2nii = pe.Node(interface=fs.preprocess.MRIConvert(), name='aparcaseg_mgz2nii')
    aparcaseg_mgz2nii.inputs.in_type = 'mgz'
    aparcaseg_mgz2nii.inputs.out_type = 'niigz'
    aparcaseg_mgz2nii.inputs.out_orientation = 'RAS'
    reorient2std = pe.Node(interface=fsl.utils.Reorient2Std(), name='reorient2std')
    diff2struct = pe.Node(interface=fsl.preprocess.FLIRT(), name='diff2struct')
    diff2struct.inputs.searchr_x = [180, 180]
    diff2struct.inputs.searchr_y = [180, 180]
    diff2struct.inputs.searchr_z = [180, 180]
    diff2struct.inputs.cost = 'mutualinfo'
    diff2struct.inputs.output_type = "NIFTI_GZ"
    convertxfm = pe.Node(interface=fsl.utils.ConvertXFM(), name='convertXFM')
    convertxfm.inputs.invert_xfm = True
    aparcaseg2diff = pe.Node(interface=fsl.preprocess.ApplyXfm(), name='aparcaseg2diff')
    aparcaseg2diff.inputs.apply_xfm = True
    aparcaseg2diff.inputs.interp = 'nearestneighbour'
    T12diff = pe.Node(interface=fsl.preprocess.ApplyXfm(), name='T12diff')
    T12diff.inputs.apply_xfm = True
    T12diff.inputs.interp = 'nearestneighbour'
    outputnode = pe.Node(niu.IdentityInterface(fields=['out_aparcaseg_diff', 'out_T1_diff']), name='outputnode')

    wf = pe.Workflow(name=name)
    wf.connect([
        (inputnode, T1_mgz2nii, [('in_T1', 'in_file')]),
        (inputnode, aparcaseg_mgz2nii, [('in_aparcaseg', 'in_file')]),
        (aparcaseg_mgz2nii, reorient2std, [('out_file', 'in_file')]),
        (inputnode, diff2struct, [('in_lowb', 'in_file')]),
        (T1_mgz2nii, diff2struct, [('out_file', 'reference')]),
        (diff2struct, convertxfm, [('out_matrix_file', 'in_file')]),
        (aparcaseg_mgz2nii, aparcaseg2diff, [('out_file', 'in_file')]),
        (inputnode, aparcaseg2diff, [('in_lowb', 'reference')]),
        (convertxfm, aparcaseg2diff, [('out_file', 'in_matrix_file')]),
        (T1_mgz2nii, T12diff, [('out_file', 'in_file')]),
        (inputnode, T12diff, [('in_lowb', 'reference')]),
        (convertxfm, T12diff, [('out_file', 'in_matrix_file')]),
        (aparcaseg2diff, outputnode, [('out_file', 'out_aparcaseg_diff')]),
        (T12diff, outputnode, [('out_file', 'out_T1_diff')])
    ])

    return wf
