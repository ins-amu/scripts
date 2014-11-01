#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include "libremesh/modelwriter.h"
#include "libremesh/defines.h"
#include "libremesh/hrtimer.h"
#include "libremesh/interface.h"

#if 1
/* Beethoven settings */
# define MESH_VERTS 500000
# define LLOYD_ITER 100
# define USE_DENSITY_FIELD true
# define DENSITY_MAX_CURVATURE 8000
# define DENSITY_CONTRAST_EXP 3.0f
# define DENSITY_SMOOTH_ITER 50
#endif

#if 0
/* Bunny settings. */
# define MESH_VERTS 20000
# define LLOYD_ITER 100
# define USE_DENSITY_FIELD true
# define DENSITY_MAX_CURVATURE 500
# define DENSITY_CONTRAST_EXP 2.0f
# define DENSITY_SMOOTH_ITER 20
#endif

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
    std::size_t density_time = 0;
    std::size_t resampling_time = 0;
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

    // Density field
    if (USE_DENSITY_FIELD)
    {
        Remesher::DensityFieldConf density_conf;
        density_conf.max_curvature = DENSITY_MAX_CURVATURE;
        density_conf.contrast_exp = DENSITY_CONTRAST_EXP;
        density_conf.smooth_iter = DENSITY_SMOOTH_ITER;
        iface.set_density_field_conf(density_conf);

        Remesher::HRTimer t_density;
        iface.exec_density_calculation();
        density_time = t_density.get_elapsed();
    }

    // Resampling
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

    // Saving
    {
        Remesher::TriangleMeshPtr mesh(iface.get_evolving_mesh());
        Remesher::ModelWriter::save_model(outfile, mesh);
    }

    std::stringstream ss;

    ss << "Cleaning took " << cleaning_time << "ms." << std::endl;
    ss << "Density took " << density_time << "ms." << std::endl;
    ss << "Resampling took " << resampling_time << "ms." << std::endl;
    ss << "Lloyd took " << lloyd_time << "ms (using "
        << (REMESHER_PARALLELIZATION * REMESHER_RELAXATION_THREADS)
        << " threads)" << std::endl;
    ss << "Total operation took " << t_total.get_elapsed()
        << "ms (" << (resampling_time + lloyd_time) << "ms "
        << "for resampling and relaxation)" << std::endl;

    std::cout << std::endl;
    std::cout << ss.str();
    std::ofstream log((outfile + ".log").c_str());
    log << ss.str();
    log.close();

    return 0;
}
