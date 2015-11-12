import unittest
import utility as ut
import numpy as np

# run tests with 
# cd scripts
# python -m unittest -v tests.test

# test subject is HCP: 151526

# class TestSurfaceBash(unittest.TestCase):
# 	def test_bash_extract_high():
# 		run extract_high lh


class TestSurfaceNodes(unittest.TestCase):
    path_s = 'tests/test_data/surface/'  # test data surface
    path_pd = 'tests/test_data/produced_data/'

    def test_fs2txt(self):
        extract_high = ut.Fs2Txt()
        extract_high.inputs.surface = path_td + 'lh.pial.asc'
        extract_high.inputs.out_file_triangles = path_pd + 'llh_triangles_high.txt'
        extract_high.inputs.out_file_vertices = path_pd + 'lh_vertices_high.txt'
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
        txt2off.inputs.out_file_off = 'tests/test_data/produced_data/lh_high.off'
        txt2off.run()
        ref_surf = open('tests/test_data/surface/lh_high.off').readlines()
        res_surf = open('tests/test_data/produced_data/lh_high.off').readlines()
        self.assertEqual(ref_surf, res_surf)

    def test_remesher(self):
        remesher = ut.Remesher()
        remesher.inputs.in_file = 'tests/test_data/surface/lh_high.off'
        remesher.inputs.out_file = 'tests/test_data/produced_data/lh_low.off'
        remesher.run()
        ref_surf = open('tests/test_data/surface/lh_low.off').readlines(2)
        res_surf = open('tests/test_data/produced_data/lh_low.off').readlines(2)
        self.assertEqual(ref_surf, res_surf)

    def test_off2txt(self):
        off2txt = ut.Off2Txt()
        off2txt.inputs.surface_off = 'tests/test_data/surface/lh_low.off'
        off2txt.inputs.out_file_vertices_txt = 'tests/test_data/produced_data/lh_vertices_low.txt'
        off2txt.inputs.out_file_triangles_txt = 'tests/test_data/produced_data/lh_triangles_low.txt'
        off2txt.run()
        ref_surf = np.loadtxt('tests/test_data/surface/lh_vertices_low.txt')
        res_surf = np.loadtxt('tests/test_data/produced_data/lh_vertices_low.txt')
        np.testing.assert_array_equal(ref_surf, res_surf)

    def test_regionmapping(self):
        rm = ut.RegionMapping()
        rm.inputs.aparc_annot = 'tests/test_data/label/lh.aparc.annot'
        rm.inputs.ref_table = 'lh_ref_table.txt'
        rm.inputs.vertices_downsampled = 'tests/test_data/surface/lh_vertices_low.txt'
        rm.inputs.triangles_downsampled = 'tests/test_data/surface/lh_triangles_low.txt'
        rm.inputs.vertices = 'tests/test_data/surface/lh_vertices_high.txt'
        rm.inputs.scripts_directory = '.'
        rm.inputs.out_file = 'tests/test_data/produced_data/region_mapping_low_not_corrected.txt'
        rm.run()
        ref_rm = np.loadtxt('tests/test_data/surface/region_mapping_low_not_corrected.txt')
        res_rm = np.loadtxt('tests/test_data/produced_data/region_mapping_low_not_corrected.txt')
        np.testing.assert_array_equal(ref_rm, res_rm)

    # def test_correctregionmapping(self):
    #     rm = ut.CorrectRegionMapping()
    #     rm.inputs.vertices = 'tests/test_data/surface/lh_vertices_low.txt'
    #     rm.inputs.triangles = 'tests/test_data/surface/lh_triangles_low.txt'
    #     rm.inputs.texture = 'tests/test_data/surface/lh_region_mapping_low_not_corrected.txt'
    #     rm.inputs.out_file_region_mapping_corr = 'tests/test_data/produced_data/region_mapping_low.txt'
    #     rm.run()
    #     ref_rm = np.loadtxt('tests/test_data/surface/region_mapping_low.txt')
    #     res_rm = np.loadtxt('tests/test_data/produced_data/region_mapping_low.txt')
    #     np.testing.assert_array_equal(ref_rm, res_rm)

    def test_checkregionmapping(self):
        rm = ut.CheckRegionMapping()
        rm.inputs.vertices = 'tests/test_data/surface/lh_vertices_low.txt'
        rm.inputs.triangles = 'tests/test_data/surface/lh_triangles_low.txt'
        rm.inputs.region_mapping = 'tests/test_data/surface/lh_region_mapping_low.txt'
        rm.inputs.scripts_directory = '.'
        rm.inputs.check = False
        rm.inputs.display= False
        rm.inputs.out_file = 'tests/test_data/produced_data/region_mapping_low.txt'
        rm.run()
        ref_rm = np.loadtxt('tests/test_data/surface/region_mapping_low.txt')
        res_rm = np.loadtxt('tests/test_data/produced_data/region_mapping_low.txt')
        np.testing.assert_array_equal(ref_rm, res_rm)

    def test_reunifybothhemisphere(self):
        rm = ut.RegionMapping()
        rm.inputs.aparc_annot = 'tests/test_data/label/lh.aparc.annot'
        rm.inputs.vertices = ['tests/test_data/surface/lh_vertices_low.txt',
                              'tests/test_data/surface/rh_vertices_low.txt']
        rm.inputs.triangles= ['tests/test_data/surface/lh_triangles_low.txt',
                              'tests/test_data/surface/rh_triangles_low.txt']
        rm.inputs.textures = ['tests/test_data/subj/lh_region_mapping_low.txt',
                              'tests/test_data/subj/rh_region_mapping_low.txt']
        rm.inputs.out_file_vertices = 'tests/test_data/produced_data/vertices.txt'
        rm.inputs.out_file_triangles = 'tests/test_data/produced_data/triangles.txt'
        rm.inputs.out_file_texture = 'tests/test_data/produced_data/region_mapping.txt'
        rm.run()
        ref_rm = np.loadtxt('tests/test_data/subj/region_mapping.txt')
        res_rm = np.loadtxt('tests/test_data/produced_data/region_mapping.txt')
        np.testing.assert_array_equal(ref_rm, res_rm)


class TestSubcorticalSurface(unittest.TestCase):
    def test_aseg2srf(self):
        pass

    def test_listsubcortical(self):
        pass

if '__name__' == '__main__':
    unittest.main()
