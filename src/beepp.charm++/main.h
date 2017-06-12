#ifndef __MAIN_H__
#define __MAIN_H__
#include <boost/shared_ptr.hpp>
#include "../bem/config_file.h"

class Main : public Chare
{

private:

    double universe_edge_length;
    double Dsolvent;
    double kappa;
    unsigned int total_num_patches;
    unsigned int total_num_charges;
    std::string output_filename;
    long start_clock;

    unsigned int default_quad_points_per_triangle;
    unsigned int default_qual_points_per_triangle;


public:

	/// Constructors ///
    Main(CkArgMsg* msg);
    Main(CkMigrateMessage* msg);

    /// Entry Methods ///
    void create_workers(CkCallback);
    void get_rhs();

    // This gets called when all goes quiet
    void quiescenceHandler();

    // gets called on every processor at start
    static void initManualLB() {
    	std::cout << "Disabling load balancing on processor " << CkMyPe() << std::endl;
    	TurnManualLBOn();
    }

    boost::shared_ptr<ConfigFile> config;

};


#define CK_TEMPLATES_ONLY
//#include "bem_fmm.def.h"
#undef CK_TEMPLATES_ONLY

#define CK_TEMPLATES_ONLY
//#include "vanilla_fmm.def.h"
#undef CK_TEMPLATES_ONLY

#endif //__MAIN_H__
