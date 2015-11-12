#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include "libremesh/modelwriter.h"
#include "libremesh/defines.h"
#include "libremesh/hrtimer.h"
#include "libremesh/interface.h"

/* Beethoven settings */
/* resampling */
# define USE_RESAMPLING false
# define MESH_VERTS 10000
/* lloyd */
# define LLOYD_ITER 20
/* density field */
# define USE_DENSITY_FIELD false
/* simplification */
# define USE_SIMPLIFICATION true
# define VERTEX_LIMIT 10000
# define FIDELITY_CHECK true
# define FIDELITY_PENALTY 10.0f
# define AREA_CHECKS true
# define AREA_PENALTY 10.0f
# define LOCAL_FEATURES true
# define LOCAL_FEATURE_ANGLE 100.0f
# define GLOBAL_FEATURE false
/* feature extraction */
# define USE_FEATURES_EXTRACTION false
# define  USE_ANGLE true
# define ANGLE 1.f
# define COMPLEX_EDGES true
# define BORDER_EDGES true;

int
main (int argc, char** argv)
{
    if (argc < 3)
    {
        std::cout << "Syntax: " << argv[0] << " <input mesh> <output mesh>"
            << std::endl;
        std::cout << "Resulting mesh is written to ./output.off" << std::endl;
        return 1;
    }
    std::string infile(argv[1]);
    std::string outfile(argv[2]);

    Remesher::HRTimer t_total;

    std::size_t cleaning_time = 0;
    std::size_t feature_time = 0;
    std::size_t resampling_time = 0;
    std::size_t simplification_time = 0;
    std::size_t lloyd_time = 0;

    // Loading, cleaning
    Remesher::Interface iface;
    {
        iface.load_model(infile);
        Remesher::HRTimer t_cleaning;
        //iface.get_reference_mesh()->scale_and_center();
        iface.clean_reference_mesh();
        iface.optimize_reference_mesh();
        cleaning_time = t_cleaning.get_elapsed();

        Remesher::TriangleMeshPtr mesh = iface.get_reference_mesh();

        std::cout << "Reference mesh has " << mesh->get_vertices().size()
            << " vertices and " << mesh->get_faces().size() << " faces"
            << std::endl;
    }

    //  Features Extraction
    if (USE_FEATURES_EXTRACTION)
    {
        Remesher::FeatureEdgesConf features_edges_conf;
        // use dihedral angle for feature extraction 
        features_edges_conf.use_angle = USE_ANGLE;
        //  angle to use
        features_edges_conf.angle = ANGLE;
        // mesh boundary
        features_edges_conf.complex_edges = COMPLEX_EDGES;
        // contour edges
        features_edges_conf.border_edges = BORDER_EDGES;
        iface.set_feature_edges_conf(features_edges_conf);

        Remesher::HRTimer t_features_edges;
        iface.exec_feature_extraction();
        feature_time = t_features_edges.get_elapsed();
    } 

    // Simplification
    if (USE_SIMPLIFICATION)
    {
        Remesher::SimplificationConf simplification_conf;
        // Vertex Limit
        simplification_conf.vertex_limit = VERTEX_LIMIT;
        // Fidelity check
        simplification_conf.perform_fidelity_checks = FIDELITY_CHECK;
        simplification_conf.fidelity_check_penalty = FIDELITY_PENALTY;
        // Vertex area checks
        simplification_conf.perform_area_checks = AREA_CHECKS;
        simplification_conf.area_check_penalty = AREA_PENALTY;
        // Check for local features
        simplification_conf.check_local_features = LOCAL_FEATURES;
        iface.set_simplification_conf(simplification_conf);
        // Check global features
        simplification_conf.keep_global_features = GLOBAL_FEATURE;

        Remesher::HRTimer t_simplification;
        iface.exec_simplification();
        simplification_time = t_simplification.get_elapsed();
    }

    // Resampling
    if (USE_RESAMPLING)
    {
        Remesher::ResamplingConf resampling_conf;
        resampling_conf.sample_amount = MESH_VERTS;
        resampling_conf.perform_decimation = true;
        iface.set_resampling_conf(resampling_conf);

        Remesher::HRTimer t_resampling;
        iface.exec_resampling();
        resampling_time = t_resampling.get_elapsed();
    }



    // Lloyd
    {
        Remesher::RelaxationConf lloyd_conf;
        lloyd_conf.iterations = LLOYD_ITER;
        iface.set_lloyd_conf(lloyd_conf);

        Remesher::HRTimer t_lloyd;
        iface.exec_lloyd();
        lloyd_time = t_lloyd.get_elapsed();
    }
    // reclean resulting mesh
    {
        iface.clean_evolving_mesh();
    }
    // Saving
    {
        Remesher::TriangleMeshPtr mesh(iface.get_evolving_mesh());
        Remesher::ModelWriter::save_model(outfile, mesh);
    }

    std::cout << "Cleaning took " << cleaning_time << "ms." << std::endl;
    std::cout << "Features extraction took " << feature_time << "ms." << std::endl;
    std::cout << "Resampling took " << resampling_time << "ms." << std::endl;
    std::cout << "Simplification took " << simplification_time << "ms." << std::endl;
    std::cout << "Lloyd took " << lloyd_time << "ms (using "
        << (REMESHER_PARALLELIZATION * REMESHER_RELAXATION_THREADS)
        << " threads)" << std::endl;
     std::cout << "Total operation took " << t_total.get_elapsed()
        << "ms (" << (resampling_time + feature_time + simplification_time + lloyd_time) << "ms "
        << "for fatures extraction and simplification and relaxation)" << std::endl;


    return 0;
}
