from nipype import SelectFiles, Node
import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import nipype.interfaces.utility as niu
import workflows as wf
import os

subject_id = 'af'
subject_directory = '/Users/timp/Work/freesurfer/'
data_directory = '/Users/timp/Work/data/'
freesurfer_directory = '/Applications/freesurfer/'
scripts_directory = '/Users/timp/Desktop/scripts/'


# inputs of the pipeline: DWI (nii, bvecs, bvals), T1, subject_id
# you can use any directory hierarchy you want for your data
# you just have to reflect it in the template
# glob syntax is allowed, for instance * ? [0-3]
sf = pe.Node(niu.IdentityInterface(
    fields=['T1', 'DWI', 'bvecs', 'bvals', 'scripts_directory', 'subject_id', 'subjects_dir', 'T1']),
             name='select_files')
# Input data in nii format:
sf.inputs.T1 = os.path.join(data_directory, subject_id, 'data/T1/T1.nii.gz')
sf.inputs.DWI = os.path.join(data_directory, subject_id, 'data/DWI/DWI.nii.gz')
sf.inputs.bvecs = os.path.join(data_directory, subject_id, 'data/DWI/bvecs')
sf.inputs.bvals = os.path.join(data_directory, subject_id, 'data/DWI/bvals')
sf.inputs.subject_id = subject_id
sf.inputs.scripts_directory = freesurfer_directory
sf.inputs.subjects_dir = subject_directory
sf.inputs.subject_id = subject_id

# if data in dcm format, please convert to nii format with mrconvert (mrtrix)
# mrconvert DWI/ dwi.nii.gz
# mrinfo DWI -export_grad_fsl bvecs bvals
# mrconvert T1/ T1.nii.gz

# instantiation pipeline
sc = wf.scripts()
# options
sc.inputs.surface.check_region_mapping.check = False
sc.inputs.surface.check_region_mapping.display = False
sc.inputs.tractography.dwi2fod.max_sh = 6
sc.inputs.tractography.tckgen.n_tracks = 1000
sc.inputs.tractography.tcksift.term_number = 500
sc.inputs.tractography.labelconfig.in_config = scripts_directory + 'fs_region.txt'
sc.inputs.tractography.labelconfig.lut_fs = freesurfer_directory + 'FreeSurferColorLUT.txt'
sc.inputs.tractography.compute_connectivity.corr_table = scripts_directory + 'correspondance_mat.txt'
sc.inputs.tractography.compute_connectivity.name_regions = scripts_directory + 'name_regions.txt'

# outputs of the pipeline: DataSink
# collecting output data: DataSink
ds = pe.Node(nio.DataSink(['lh', 'rh']), name='sinker')
# Directory where data are going to be outputed
ds.inputs.base_directory = os.path.join('/Users/timp/Work/processed_data', subject_id)

# Main workflow
main_wf = pe.Workflow(name='main')
main_wf.connect([
    (sf, sc, [('T1', 'inputnode.T1'),
              ('DWI', 'inputnode.DWI'),
              ('bvecs', 'inputnode.bvecs'),
              ('bvals', 'inputnode.bvals'),
              ('scripts_directory', 'inputnode.scripts_directory'),
              ('subjects_dir', 'inputnode.subjects_dir'),
              ('subject_id', 'inputnode.subject_id')]),
    (sc, ds, [('outputnode.weights', 'connectivity.@weights'),  # i.e. in connectivity directory
              ('outputnode.tract_lengths', 'connectivity.@tract_lengths'),
              ('outputnode.areas', 'connectivity.@areas'),
              ('outputnode.centres', 'connectivity.@centres'),
              ('outputnode.average_orientations', 'connectivity.@average_orientations'),
              ('outputnode.texture', 'surface.@region_mapping'),
              ('outputnode.vertices', 'surface.@vertice'),
              ('outputnode.triangles', 'surface.@triangles')])
])

# Directory for temporary data
main_wf.base_dir = '/Users/timp/Work/tmp_pip/'
# Output a graph of the processing steps
main_wf.write_graph()
# Run the pipeline
main_wf.run()
