from nipype import SelectFiles, Node
import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import workflows as wf
import os

subject_id = 'af'

# inputs of the pipeline: DWI (nii, bvecs, bvals), T1, subject_id
# you can use any directory hierarchy you want for your data
# you just have to reflect it in the template
# glob syntax is allowed, for instance * ? [0-3]

# Input data in nii format:
templates = dict(T1="{subject_id}/data/T1/*.nii.gz",
                 DWI="{subject_id}/data/DWI/*.nii*",
                 bvecs="{subject_id}/data/DWI/bvecs",
                 bvals="{subject_id}/data/DWI/bvals",
                 subject_id="{subject_id}",
                 scripts_directory=os.getcwd(),
                 subjects_dir="/Users/timp/Work/freesurfer/")
sf = Node(SelectFiles(templates), 'select_files')
sf.inputs.subject_id = subject_id
# Base directory for where are the data
sf.inputs.base_directory = '/Users/timp/Work/data/'

# if data in dcm format, please convert to nii format with mrconvert (mrtrix)
# mrconvert DWI/ dwi.nii.gz
# mrinfo DWI -export_grad_fsl bvecs bvals
# mrconvert T1/ T1.nii.gz

# instantiation pipeline
sc = wf.scripts()
# options
sc.inputs.surface.check_region_mapping.check = False
sc.inputs.surface.check_region_mapping.display = False
sc.inputs.tractography.dwi2fod.max_sh = 8
sc.inputs.tractography.tckgen.n_tracks = 1000
sc.inputs.tractography.tcksift.term_number = 500

# outputs of the pipeline: DataSink
# collecting output data: DataSink
ds = pe.Node(nio.DataSink(), name='sinker')
# Directory where data are going to be outputed
ds.inputs.base_directory = '/Users/timp/Work/processed_data/{subject_id}/'



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

