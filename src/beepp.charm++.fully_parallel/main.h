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
	std::string output_filename;
    long start_clock;


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
//#include "parallel_fmm_octree.def.h"
#undef CK_TEMPLATES_ONLY

#define CK_TEMPLATES_ONLY
//#include "vanilla_fmm_evals.def.h"
#undef CK_TEMPLATES_ONLY

#define CK_TEMPLATES_ONLY
//#include "vanilla_fmm_worker.def.h"
#undef CK_TEMPLATES_ONLY

#endif //__MAIN_H__
