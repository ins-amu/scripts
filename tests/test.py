import unittest
import utility as ut
import numpy as np

# run tests with 
# cd scripts
# python -m unittest -v tests.test

#class TestSurfaceBash(unittest.TestCase):
#	def test_bash_extract_high():
#		run extract_high lh

class TestSurfaceNodes(unittest.TestCase):
	def test_fs2txt(self):
		extract_high = ut.Fs2Txt()
		extract_high.inputs.surface = 'tests/test_data/surface/lh.pial.asc'
		extract_high.inputs.out_file_triangles = 'tests/test_data/produced_data/lh_triangles_high.txt'
		extract_high.inputs.out_file_vertices = 'tests/test_data/produced_data/lh_vertices_high.txt'
		extract_high.run()
		ref_tri = np.loadtxt('tests/test_data/surface/lh_triangles_high.txt')
		ref_vert = np.loadtxt('tests/test_data/surface/lh_vertices_high.txt')
		res_tri = np.loadtxt('tests/test_data/produced_data/lh_triangles_high.txt')
		res_vert = np.loadtxt('tests/test_data/produced_data/lh_vertices_high.txt')
		np.testing.assert_array_equal(ref_tri, res_tri)
		np.testing.assert_array_equal(ref_vert, res_vert)

	def test_txt2off(self):
		txt2off = ut.Txt2Off()
		txt2off.inputs.vertices = 'tests/test_data/surface/lh_vertices_high.txt'
		txt2off.inputs.triangles = 'tests/test_data/surface/lh_triangles_high.txt'
		txt2off.out_file_off = 'tests/test_data/produced_data/lh_high.off'
		txt2off.run()
		ref_surf = np.loadtxt('tests/test_data/surface/lh_high.off', skiprows=1)
		res_surf = np.loadtxt('tests/test_data/produced_data/lh_high.off', skiprows=1)
		np.testing.assert_array_equal(ref_surf, res_surf)

	def test_regionmapping(self):
		rm = ut.RegionMapping()
		rm.inputs.aparc_annot = 'tests/test_data/label/lh.aparc.annot'
		rm.inputs.ref_table = 'lh_ref_table.txt'
		rm.inputs.vertices_downsampled = 'tests/test_data/surface/lh_vertices_low.txt'
		rm.inputs.triangles_downsampled = 'tests/test_data/surface/lh_triangles_low.txt'
		rm.inputs.vertices = 'tests/test_data/surface/lh_vertices_high.txt'
		rm.inputs.scripts_directory = '.'
		rm.inputs.out_file = 'tests/test_data/produced_data/region_mapping.txt'
		rm.run()
		ref_rm = np.loadtxt('tests/test_data/subj/region_mapping.txt')
		res_rm = np.loadtxt('tests/test_data/produced_data/region_mapping.txt')
		np.testing.assert_array_equal(ref_rm, res_rm)

if '__name__'=='__main__':
	unittest.main()
