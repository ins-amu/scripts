import unittest
import utility as ut
import numpy as np
import workflows as wf


# run tests with:
# cd scripts
# python -m unittest -v tests.test


class TestSurface(unittest.TestCase):
    """
    Test for surface workflow

    Test data are from HCP subject 151526
    """
    path_s = 'tests/test_data/surface/'  # test data surface
    path_pd = 'tests/test_data/produced_data/'  # test produced data
    path_d = 'tests/test_data/subj/'  # test data subj

    def test_fs2txt(self):
        extract_high = ut.Fs2Txt()
        extract_high.inputs.surface = self.path_s + 'lh.pial.asc'
        extract_high.inputs.out_file_triangles = self.path_pd + 'lh_triangles_high.txt'
        extract_high.inputs.out_file_vertices = self.path_pd + 'lh_vertices_high.txt'
        extract_high.run()
        ref_tri = np.loadtxt(self.path_s + 'lh_triangles_high.txt')
        ref_vert = np.loadtxt(self.path_s + 'lh_vertices_high.txt')
        res_tri = np.loadtxt(self.path_pd + 'lh_triangles_high.txt')
        res_vert = np.loadtxt(self.path_pd + 'lh_vertices_high.txt')
        np.testing.assert_array_equal(ref_tri, res_tri)
        np.testing.assert_array_equal(ref_vert, res_vert)

    def test_txt2off(self):
        txt2off = ut.Txt2Off()
        txt2off.inputs.vertices = self.path_s + 'lh_vertices_high.txt'
        txt2off.inputs.triangles = self.path_s + 'lh_triangles_high.txt'
        txt2off.inputs.out_file_off = self.path_pd + 'lh_high.off'
        txt2off.run()
        ref_surf = open(self.path_s + 'lh_high.off').readlines()
        res_surf = open(self.path_pd + 'lh_high.off').readlines()
        self.assertEqual(ref_surf, res_surf)

    def test_remesher(self):
        remesher = ut.Remesher()
        remesher.inputs.in_file = self.path_s + 'lh_high.off'
        remesher.inputs.out_file = self.path_pd + 'lh_low.off'
        remesher.run()
        f = open(self.path_s + 'lh_low.off')
        f.readline()
        ref_surf = f.readline()
        g = open(self.path_pd + 'lh_low.off')
        g.readline()
        res_surf = g.readline()
        self.assertEqual(ref_surf, res_surf)

    def test_off2txt(self):
        off2txt = ut.Off2Txt()
        off2txt.inputs.surface_off = self.path_s + 'lh_low.off'
        off2txt.inputs.out_file_vertices_txt = self.path_pd + 'lh_vertices_low.txt'
        off2txt.inputs.out_file_triangles_txt = self.path_pd + 'lh_triangles_low.txt'
        off2txt.run()
        ref_surf = np.loadtxt(self.path_s + 'lh_vertices_low.txt')
        res_surf = np.loadtxt(self.path_pd + 'lh_vertices_low.txt')
        np.testing.assert_array_equal(ref_surf, res_surf)

    def test_regionmapping(self):
        rm = ut.RegionMapping()
        rm.inputs.vertices_downsampled = self.path_s + 'lh_vertices_low.txt'
        rm.inputs.vertices = self.path_s + 'lh_vertices_high.txt'
        rm.inputs.aparc_annot = 'tests/test_data/label/lh.aparc.annot'
        rm.inputs.ref_table = 'lh_ref_table.txt'
        rm.inputs.scripts_directory = '.'
        rm.inputs.out_file = self.path_pd + 'lh_region_mapping_low_not_corrected.txt'
        rm.run()
        ref_rm = np.loadtxt(self.path_s + 'lh_region_mapping_low_not_corrected.txt')
        res_rm = np.loadtxt(self.path_pd + 'lh_region_mapping_low_not_corrected.txt')
        np.testing.assert_array_equal(ref_rm, res_rm)


    def test_correctandcheckregionmapping(self):
        rm = ut.CorrectRegionMapping()
        rm.inputs.vertices = self.path_s + 'lh_vertices_low.txt'
        rm.inputs.triangles = self.path_s + 'lh_triangles_low.txt'
        rm.inputs.texture = self.path_s + 'lh_region_mapping_low_not_corrected.txt'
        rm.inputs.out_file_region_mapping_corr = self.path_pd + '/region_mapping_low.txt'
        rm.run()
        rm = ut.CheckRegionMapping()
        rm.inputs.vertices = self.path_s + 'lh_vertices_low.txt'
        rm.inputs.triangles = self.path_s + 'lh_triangles_low.txt'
        rm.inputs.region_mapping = self.path_pd + 'region_mapping_low.txt'
        rm.inputs.scripts_directory = '.'
        rm.inputs.check = False
        rm.inputs.display = False
        rm.inputs.out_file = self.path_pd + 'region_mapping_low.txt'
        rm.run()
        ref_rm = np.loadtxt(self.path_s + 'lh_region_mapping_low.txt')
        res_rm = np.loadtxt(self.path_pd + 'region_mapping_low.txt')
        np.testing.assert_array_equal(ref_rm, res_rm)

    def test_reunifybothhemisphere(self):
        rm = ut.ReunifyBothHemisphere()
        rm.inputs.vertices = [self.path_s + 'lh_vertices_low.txt',
                              self.path_s + 'rh_vertices_low.txt']
        rm.inputs.triangles = [self.path_s + 'lh_triangles_low.txt',
                               self.path_s + 'rh_triangles_low.txt']
        rm.inputs.textures = [self.path_s + 'lh_region_mapping_low.txt',
                              self.path_s + 'rh_region_mapping_low.txt']
        rm.inputs.out_file_vertices = self.path_pd + 'vertices.txt'
        rm.inputs.out_file_triangles = self.path_pd + 'triangles.txt'
        rm.inputs.out_file_texture = self.path_pd + 'region_mapping.txt'
        rm.run()
        ref_rm = np.loadtxt(self.path_d + 'region_mapping.txt')
        res_rm = np.loadtxt(self.path_pd + 'region_mapping.txt')
        np.testing.assert_array_equal(ref_rm, res_rm)
        ref_verts = np.loadtxt(self.path_d + 'vertices.txt')
        res_verts = np.loadtxt(self.path_pd + 'vertices.txt')
        np.testing.assert_array_equal(ref_verts, res_verts)
        ref_tri = np.loadtxt(self.path_d + 'triangles.txt')
        res_tri = np.loadtxt(self.path_pd + 'triangles.txt')
        np.testing.assert_array_equal(ref_tri, res_tri)


class TestSubcorticalSurface(unittest.TestCase):
    def test_aseg2srf(self):
        pass

    def test_listsubcortical(self):
        pass


class TestPreprocess(unittest.TestCase):
    """
    Test data for preprocessing workflow

    Test data are from subject af
    """

    def test_mrconvert(self):
        pass


class TestConnectivity(unittest.TestCase):
    """
    Test data for preprocessing workflow

    Test data are from subject af
    """

    # use mrtrix test data
    def test_dwi_extract_lowb(self):
        pass

    def test_fivett2gmwmi(self):
        pass

    def test_compute_connectivity(self):
        pass


class TestWorkflows(unittest.TestCase):
    """
    Test instantiation of the workflows
    """

    def test_scripts(self):
        wf.scripts()

    def test_coregistration(self):
        wf.coregistration()

    def test_tractography(self):
        wf.tractography()

    def test_preprocess(self):
        wf.preprocess()

    def test_subcortical(self):
        wf.subcorticalsurface()

    def test_surface(self):
        wf.surface()


class TestCoregistration(unittest.TestCase):
    # TODO finish implementation
    def test_T1_mgz2nii(self):
        pass


class TestSurfaceCorrectness(unittest.TestCase):
    def test_single_component(self):
        verts = np.loadtxt('tests/test_data/surface/lh_vertices_low.txt')
        tri = np.loadtxt('tests/test_data/surface/lh_triangles_low.txt').astype(int)
        rm = np.loadtxt('tests/test_data/surface/lh_region_mapping_low.txt')
        crm = ut.CheckRegionMapping()
        res = crm.calculate_connected(rm, verts, tri)
        tt = np.sum(np.array(res)[:, 1])
        self.assertEqual(tt, 0)
        return res


    def test_surface_downsampling_area(self):
        verts = np.loadtxt('tests/test_data/surface/lh_vertices_low.txt')
        tri = np.loadtxt('tests/test_data/surface/lh_triangles_low.txt').astype(int)
        tri_area_low = ut.ComputeConnectivityFiles()
        res_low = np.sum(tri_area_low.compute_triangle_areas(verts, tri))
        verts_high = np.loadtxt('tests/test_data/surface/lh_vertices_high.txt')
        tri_high = np.loadtxt('tests/test_data/surface/lh_triangles_high.txt').astype(int)
        tri_area_high = ut.ComputeConnectivityFiles()
        res_high = np.sum(tri_area_high.compute_triangle_areas(verts_high, tri_high))
        self.assertGreater(res_low, .8*res_high)


    def test_surface_sphere(self):
        pass


if '__name__' == '__main__':
    unittest.main()
