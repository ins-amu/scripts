import nipype.interfaces.fsl as fsl
import nipype.interfaces.mrtrix as mrt
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as niu
import nipype.pipeline.engine as pe
from nipype.workflows.dmri.fsl.artifacts import ecc_pipeline 
import scripts.utility as su 

def ReconAll():


def Surface(name='surface'):
    """ Surface workflow
    Generate vertices.txt, triangles.txt and region_mapping.txt
    """

    inputnode = pe.Node(interface=niu.IdentityInterface(fields=['pial_rh', 'annot_rh', 'ref_tables_rh', 'pial_lh', 'annot_lh', 'ref_tables_lh','rh', 'lh']), name='inputnode')
    
    inputnode.inputs.rh ='rh'
    inputnode.inputs.lh ='lh'
    inputnode.inputs.pial_rh = "/disk3/Work/Processed_data/freesurfer/100408/surf/rh.pial"
    inputnode.inputs.annot_rh = "/disk3/Work/Processed_data/freesurfer/100408/label/rh.aparc.annot"
    inputnode.inputs.ref_tables_rh = "/home/tim/Work/Models/processing/scripts_nathan/rh_ref_table.txt"
    inputnode.inputs.pial_lh = "/disk3/Work/Processed_data/freesurfer/100408/surf/lh.pial"
    inputnode.inputs.annot_lh = "/disk3/Work/Processed_data/freesurfer/100408/label/lh.aparc.annot"
    inputnode.inputs.ref_tables_lh = "/home/tim/Work/Models/processing/scripts_nathan/lh_ref_table.txt"


    pial2asc = pe.MapNode(interface=su.MRIsConvert(), name='pial2asc', )
    pial2asc.inputs.out_datatype ='asc'
    pial2asc.inputs.normal = True
    extract_high = pe.Node(interface=niu.Function(input_names=['surface', 'rl'],
                                                  output_names=['vertices_high', 'triangles_high'],
                                                  function=su.extract_high), name='extract_high')
    txt2off = pe.Node(interface=niu.Function(input_names=['vertices', 'triangles', 'rl'],
                                             output_names=['out_file'],
                                             function=su.txt2off),name='txt2off')
    remesher = pe.Node(interface=su.Remesher(), name='remesher')
    off2txt = pe.Node(interface=niu.Function(input_names=['surface', 'rl'],
                                             output_names=['vertices_low', 'triangles_low'],
                                             function=su.off2txt), name='off2txt')
    region_mapping = pe.Node(interface=su.RegionMapping(),name='region_mapping')
    correct_region_mapping = pe.Node(interface=niu.Function(input_names=['region_mapping_not_corrected', 'vertices', 'triangles', 'rl', 'region_mapping_corr'], 
                                                            output_names = ['region_mapping_low'],
                                                            function=su.correct_region_mapping),name='correct_region_mapping')
    check_region_mapping = pe.Node(interface=su.CheckRegionMapping(), name='check_region_mapping')
    reunify_both_regions = pe.Node(interface=niu.Function(input_names = ['rh_region_mapping', 'lh_region_mapping', 'rh_vertices', 'lh_vertices', 'rh_triangles', 'lh_triangles'],
                                                          output_names = ['out_files'],
                                                          function = su.reunify_both_regions), name='reunify_both_regions')


    wfrh = pe.Workflow(name='wfrh')
    wfrh.connect([
        (pial2asc, extract_high, [('converted','surface')]),
        (extract_high, txt2off, [('vertices_high','vertices'),
                                 ('triangles_high','triangles')]),
        (extract_high, region_mapping,[('vertices_high','vertices_high')]),
        (txt2off, remesher, [('out_file','in_file')]),
        (remesher, off2txt, [('out_file','surface')]),
        (off2txt, region_mapping, [('vertices_low','vertices_low'),
                                   ('triangles_low','triangles_low')]),
        (off2txt, correct_region_mapping, [('vertices_low','vertices'),
                                           ('triangles_low','triangles')]),
        (off2txt, check_region_mapping, [('vertices_low','vertices_low'),
                                         ('triangles_low','triangles_low')]),
        (region_mapping, correct_region_mapping,[('out_file','region_mapping_not_corrected')]),
        (correct_region_mapping, check_region_mapping, [('region_mapping_low','region_mapping_low')])
        ])

    wflh = wfrh.clone(name='wflh')

    wf = pe.Workflow(name='wf_surf')
    wf.connect([
        (inputnode, wfrh,[('pial_rh','pial2asc.in_file')]),
        (inputnode, wflh,[('pial_lh','pial2asc.in_file')]),
        (inputnode, wfrh, [('annot_rh', 'region_mapping.aparc_annot'),
                       ('ref_tables_rh','region_mapping.ref_tables'),
                       ('rh','region_mapping.rl')]),
        (inputnode, wflh, [('annot_lh', 'region_mapping.aparc_annot'),
                       ('ref_tables_lh','region_mapping.ref_tables'),
                       ('lh','region_mapping.rl')]),
        (inputnode, wfrh, [('rh','correct_region_mapping.rl')]),
        (inputnode, wflh, [('lh','correct_region_mapping.rl')]),
        (inputnode, wfrh, [('rh','extract_high.rl')]),
        (inputnode, wflh, [('lh','extract_high.rl')]),
        (inputnode, wfrh, [('rh','txt2off.rl')]),
        (inputnode, wflh, [('lh','txt2off.rl')]),
        (inputnode, wfrh,[('rh','off2txt.rl')]),
        (inputnode, wflh,[('lh','off2txt.rl')]),
        (wfrh,reunify_both_regions,[('check_region_mapping.region_mapping_low', 'rh_region_mapping')]),
        (wflh,reunify_both_regions,[('check_region_mapping.region_mapping_low', 'lh_region_mapping')]),
        (wfrh,reunify_both_regions,[('off2txt.vertices_low', 'rh_vertices')]),
        (wflh,reunify_both_regions,[('off2txt.vertices_low', 'lh_vertices')]),
        (wfrh,reunify_both_regions,[('off2txt.triangles_low', 'rh_triangles')]),
        (wflh,reunify_both_regions,[('off2txt.triangles_low', 'lh_triangles')])
        ])

    wf.run()
    
    return wf
    inputnode = pe.Node(interface=niu.IdentityInterface(fields=['pial']), name='inputnode')
    pial2asc = pe.MapNode(interface=fs.utils.MRIsConvert(), name='pial2asc',
            iterfield = ['surface', 'rl'])
    pial2asc.
    extract_high = pe.Node(interface=niu.Function(input_names=['surface', 'rl'],
                                                  output_names=['vertices_high', 'triangles_high'],
                                                  function=su.extract_high)
    txt2off = pe.Node(interface=niu.Function(input_names=['vertices', 'triangles', 'rl'],
                                             output_names=['high.off'],
                                             function=su.txt2off)
    remesher = pe.Node(interface=niu.Remesher(), name='remesher') 
    off2txt = pe.Node(interface=niu.Function(input_names=['surface', 'rl'],
                                             output_names=['vertices_low', 'triangles_low'],
                                             function=su.off2txt)
    region_mapping = pe.
    correct_region_mapping = pe.Node(interface=niu.Function(input_names=[
        'region_mapping_not_corrected', 'vertices', 'triangles', 'rl', 'region_mapping_corr'), 
        output_names = ['region_mapping_low'],
        function=su.correct_region_mapping)
    check_region_mapping = pe.Node(interface=niu.CheckRegionMapping(), name='check_region_mapping') 
    reunify_both_regions




def SubcorticalSurface(name="subcorticalsurfaces"):
    """ extraction of the subcortical surfaces from FreeSurfer"""
    inputnode = pe.Node(interface=niu.IdentityInterface(fields=['in_subject_id']), name='inputnode')
    aseg2srf = pe.Node(interface=su.Aseg2Srf(), name='aseg2srf')
    list_subcortical = pe.MapNode(interface=su.ListSubcortical(), name=list_subcortical, iterfield=['in_file'])
    outputnode = pe.MapNode(interface=niu.IdentityInterface(fields=['out_fields']), name='outputnode')

    wf = pr.Workflow(name=name)
    wf.connect([
        (inputnode, aseg2srf, [('in_subject_id', 'in_subject_id')]),
        (aseg2srf, list_subcortical, [('out_files', 'in_file')]),
        (list_subcortical, outputnode, [('out_files', 'in_file')])])

    return wf





def Connectivity(name="connectivity"):
    inputnode = pe.Node(interface=niu.IdentityInterface(fields=['in_file']), name='inputnode')
    convert_dicom2nii = pe.Node(interface=mrt.MRConvert(), name='convert_dicom2nii')
    extract_bvecs_bvals = pe.Node(interface=mtr3.utils.MRInfo(), name='extract_bvecs_bvals')
# ecc = ecc_pipeline()
    convert_nii2dwi = pe.Node(interface=mrt.MRConvert(), name='convert_nii2dwi')
    create_mask = pe.Node(interface=mrt3.Dwi2Mask(), name='create_mask')
    dwi_extract_lowb = .Node(interface=mrt3.DwiExtract(), name='dwi_extract_lowb')
    lowb_mif2lowb_nii = pe.Node(interface=mrt.MRConvert(), name='lowb_mif2lowb_nii')
    cor = Coregistration()
    dwi2response = pe.Node(interface=mrt3.preprocess().Dwi2Response, name='dwi2response')
    dwi2fod = pe.Node(interface=mrt3.preprocess().Dwi2Fod, name='dwi2fod')
    dwi2fod.inputs.lmax = 8
    act_anat_prepare_fsl = pe.node(interface=mrt3.ActAnatPrepareFSL(), name='act_anat_prepare_fsl')
    5tt2gmwmi = pe.Node(interface=mrt3.5tt2Gmwmi(), name='5tt2gmwmi')
    tckgen = pe.node(interface=mrt3.tracking.Tckgen(), name='tckgen')
    tckgen.inputs.unidirectional = True
    tckgen.inputs.seed_gmwmi = 'iFOD2'
    tckgen.inputs.maxlength = 150
    tckgen.inputs.num = 5000000
    tcksift = pe.Node(interface=mrt3.TckSift(), neame='tcksift')
    tcksift.inputs.term_number = 2500000
    labelconfig = pe.Node(interface=mrt3.utils.LabelConfig(), name='labelconfig')
    tck2connectome_weights = pe.Node(interface=mrt3.tracking.Tck2Connectome(), name='tck2connectome_weights')
    tck2connectome_weights.inputs.assignement_radial_search = 2
    tck2connectome_tract_lengths = pe.Node(interface=mrt3.tracking.Tck2Connectome(), name='tck2connectome_tract_lengths')
    tck2connectome_tract_lengths.inputs.assignement_radial_search = 2
    tck2connectome_tract_lengths.inputs.zero_diagonal = True
    tck2connectome_tract_lengths.inputs.metric = 'meanlength'
    compute_connectivity = pe.Node(interface=su.ComputeConnectivityFiles(), name='compute_connecitivty')
    outputnode = pe.Node(interface=nit.IdentityInterface(fields=('out_weights', 'out_tract_lengths']), name='outputnode')

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
        (inputnode, cor, [('in_T1.mgz', 'in_T1.mgz')]),
        (inputnode, cor, [('in_aparcaseg.mgz', 'in_aparcaseg.mgz')]),
        (lowb_mif2lowb_nii, cor, [('out_file', 'in_aparcaseg.mgz')]),
        (convert_nii2dwi, dwi2response, [('out_file', 'dwi')]),
        (create_mask, dwi2response, [('out_file', 'mask')]),
        (convert_nii2dwi, dwi2fod, [('out_file', 'dwi')]),
        (dwi2response, dwi2fod, [('response', 'response')]),
        (create_mask, dwi2fod, [('out_file', 'mask')]),
        (cor, act_anat_prepare_fsl, [('out_T1_diff', 'in_file')]),
        (act_anat_prepare_fsl, 5tt2gmwmi, [('out_file', 'in_file')]),
        (dwi2fod, tckgen, [('sh_out_file', 'source')]),
        (5tt2gmwmi, tckgen, [('out_file', 'seed_gmwmi')]),
        (tckgen, tcksift, [('tracks', 'tracks')]),
        (dwi2fod, tcksift, [('sh_out_file', 'source')]),
        (cor, labelconfig, [('out_aparcaseg_diff', 'labels_in')]),
        (inputnode, labelconfig, [('fs_region.txt', 'config')]),
        (inputnode, labelconfig, [('FreeSurferColorLUT.txt', 'lut_freesurfer')]),
        (tcksift, tck2connectome_weights, [('out_file', 'tracks_in')]), 
        (labelconfig, tck2connectome_weights, [('labels_out', 'nodes_in')]),
        (tcksift, tck2connectome_tract_lengths, [('out_file', 'tracks_in')]), 
        (labelconfig, tck2connectome_tract_lengths, [('labels_out', 'nodes_in')]),
        (tck2connectome_weights, compute_connectivity, [('connectome_out', 'in_file')]),
        (tck2connectome_tract_lengths, compute_connectivity, [('connectome_out', 'in_file')]),
        (compute_connectivity, output_node, [('weights.txt', 'out_weights.txt'),
                                             ('tract_lengths.txt', 'out_tract_lengths.txt')])])


    return wf


        
        


def Coregistration(name='coregistration'):
    inputnode = pe.Node(niu.IdentityInterface(fields=['in_T1.mgz', 'in_lowb.nii',
                        'in_aparcaseg.mgz']), name='inputnode')
    T1_mgz2nii = pe.Node(interface=fs.preprocess.MRIConvert(), name='T1_mgz2nii')
    T1_mgz2nii.inputs.in_type = 'mgz'
    T1_mgz2nii.inputs.out_type = 'nii'
    T1_mgz2nii.inputs.out_orientation = 'RAS'
    aparcaseg_mgz2nii = pe.Node(interface=fs.preprocess.MRIConvert(), name='aparcaseg_mgz2nii')
    aparcaseg_mgz2nii.inputs.in_type = 'mgz'
    aparcaseg_mgz2nii.inputs.out_type = 'nii'
    aparcaseg_mgz2nii.inputs.out_orientation = 'RAS'
    reorient2std = pe.Node(interface=fsl.utils.Reorient2Std(), name='reorient2std')
    diff2struct = pe.Node(interface=fsl.preprocess.FLIRT(), name='diff2struct')
    diff2struct.inputs.searchr_x = [180, 180]
    diff2struct.inputs.searchr_y = [180, 180]
    diff2struct.inputs.searchr_z = [180, 180]
    diff2struct.inputs.cost = 'mutualinfo'
    convertxfm = pe.Node(interface=fsl.utils.ConvertXFM(), name='convertXFM')
    convertxfm.inputs.invert_xfm = True
    aparcaseg2diff = pe.Node(interface=fsl.preprocess.ApplyXFM(), name='aparcaseg2diff')
    aparcaseg2diff.inputs.apply_xfm = True
    aparcaseg2diff.inputs.interp = 'nearestneighbour'
    T12diff = pe.Node(interface=fsl.preprocess.ApplyXFM(), name='T12diff')
    T12diff.inputs.apply_xfm = True
    T12diff.inputs.interp = 'nearestneighbour'
    outputnode = pe.Node(niu.IdentityInterface(fields=['out_aparcaseg_diff', 'out_T1_diff'], name='outputnode'))
    
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







