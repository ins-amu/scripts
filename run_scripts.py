from nipype import SelectFiles, Node
import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import workflows as wf
import os

# inputs of the pipeline: DWI, T1, subject_id
# you can use any directory hierarchy you want for your data
# you just have to reflect it in the templace
# glob syntax is allowed, for instance * ? [0-3]
subject_id = 'af'

# Example dcm:
templates = dict(T1="{subject_id}/data/T1/*.nii.gz",
                 DWI="{subject_id}/data/DWI/*.nii*",
                 bvecs="{subject_id}/data/DWI/bvecs",
                 bvals="{subject_id}/data/DWI/bvals",
                 subject_id="{subject_id}",
                 scripts_directory=os.getcwd(),
                 subjects_dir="/Users/timp/Work/freesurfer/")
# Example nii:
# templates = dict(T1="{subject_id}/T1/*.nii",
#                  DWI="{subject_id}/DWI/*.nii",
#                  subject_id="{subject_id}")

sf = Node(SelectFiles(templates), 'select_files')
sf.inputs.subject_id = subject_id
sf.inputs.base_directory = '/Users/timp/Work/data/'

# options for the pipeline

# outputs of the pipeline: DataSink

ds = pe.Node(nio.DataSink(), name='sinker')
ds.inputs.base_directory = '/Users/timp/Work/processed_data/{subject_id}/'

sc = wf.scripts()
main_wf = pe.Workflow(name='main')
main_wf.connect([
    (sf, sc, [('T1', 'inputnode.T1'),
              ('DWI', 'inputnode.DWI'),
              ('bvecs', 'inputnode.bvecs'),
              ('bvals', 'inputnode.bvals'),
              ('scripts_directory', 'inputnode.scripts_directory'),
              ('subjects_dir', 'inputnode.subjects_dir'),
              ('subject_id', 'inputnode.subject_id')]),
    (sc, ds, [('outputnode.weights', 'connectivity.@weights'),
              ('outputnode.tract_lengths', 'connectivity.@tract_lengths'),
              ('outputnode.areas', 'connectivity.@areas'),
              ('outputnode.centres', 'connectivity.@centres'),
              ('outputnode.average_orientations', 'connectivity.@average_orientations'),
              ('outputnode.texture', 'surface.@region_mapping'),
              ('outputnode.vertices', 'surface.@vertice'),
              ('outputnode.triangles', 'surface.@triangles')])
])

main_wf.base_dir = '/Users/timp/Work/tmp_pip/'
main_wf.write_graph()
main_wf.run()
