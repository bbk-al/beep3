// The main function for BEEP
#include "beep.h"
#include <vector>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
    long start = myclock();
    // Usage: bem_solver config.xml

    // get config filename as a string
    std::string config_filename;
    int cmdline_quad_points = -1;
    int cmdline_qual_points = -1;
    int cmdline_nbsize = -1;
    double cmdline_kappa = -1;
    bool use_planar;


    // positional options: input config (xml) filename is one of these
    po::positional_options_description p;
    p.add("input-file", -1);

    // Declare the supported options.
    po::options_description generic("Allowed options");
    generic.add_options()
        ("help", "Display help") 
        ("qual", po::value<int>(&cmdline_qual_points), "set qualocation rule (num pts per triangle on source patches: 0,1,4,7)")
        ("quad", po::value<int>(&cmdline_quad_points), "set quadrature rule (num pts per triangle on target patches: 0,1,4,7)")
        ("nbsize", po::value<int>(&cmdline_nbsize), "set BEM neighbourhood size")
        ("kappa", po::value<double>(&cmdline_kappa), "set BEM neighbourhood size")
        ("planar", po::value<bool>(&use_planar)->default_value(false), "force use of planar triangles")
    ;

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("input-file", po::value<std::string>(&config_filename)->default_value(""), "xml format configuration input file")
        ;
        
    po::options_description opts;
    opts.add(generic).add(hidden);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
            options(opts).positional(p).run(), vm);
    po::notify(vm);

    // TODO: put something useful here
    if (vm.count("help")) {
        std::cout << generic << "\n";
        return 0;
    }
    
    // handle no config file
    if (config_filename == "")
    {
        // No input file
        std::cerr << "No input file specified." << std::endl;
        return 1;
    }
    
    try
    {
        // load config file
        ConfigFile config(config_filename);

        // Create a bem-solver, using the config
        BEEP BEM_solver(config, cmdline_quad_points, cmdline_qual_points, cmdline_nbsize, cmdline_kappa, use_planar);

        // run BEM electrostatics -- 1e-6 gmres tolerance, max 100 iterations
        BEM_solver.solve(1e-6, 100);
        
        BEM_solver.calculate_energies();
        BEM_solver.calculate_forces();

        // write the fh values to output file
        BEM_solver.write_fh(config.output_file);

    }
    catch (ConfigFile::ParseException)
    {
        std::cerr << "Failed to parse config file." << std::endl;
        return 1;
    }
    catch (Mesh::MeshError)
    {
        std::cerr << "Failed to load a Mesh." << std::endl;
        return 1;
    }

    std::cout << "Total time: " << (myclock() - start)/1000. << std::endl;

    return 0;
}
