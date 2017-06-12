// things that appear in the .decl files and therefore need to be included first
#include "vanilla_fmm_prerequisites.h"

// Other Chare types
#include "main.decl.h"
#include "vanilla_fmm.decl.h"
#include "vanilla_fmm.h"
#include "main.decl.h"
#include "main.h"
#include "../common/useful_clock.h"

// 'Normal' C++ headers
#include <string>
#include <boost/lexical_cast.hpp>
#include "../fmm/eval_pt.h"

#define KAPPA 0.1

// FMM Max objects per node
#ifdef OPENCL
static const unsigned int default_max_objects = 1000;
#else
static const unsigned int default_max_objects = 1000;
#endif

/* readonly */ CProxy_Main MainProxy;
/* readonly */ CProxy_VanillaFMM<18,18,300> VanillaFMMProxy;
/* readonly */ CProxy_OpenCL_NodeGroup OpenCL_NodeGroupProxy;

typedef std::numeric_limits< double > dbl;

int get_random_number(int low, int high)
{
    int range=(high-low)+1;
    return low+int(range*double(rand())/(RAND_MAX + 1.0));
}

int random_sign()
{
    return (rand() >= RAND_MAX/2) ? -1 : 1;
}

// Entry point of Charm++ application
Main::Main(CkArgMsg* msg)
{
    start_time = myclock();
    CkPrintf("Running Parallel FMM using %d processors.\n", CkNumPes());

    // Set the mainProxy readonly to point to a
    // proxy for the Main chare object (this
    // chare object).
    MainProxy = thisProxy;
    OpenCL_NodeGroupProxy = CProxy_OpenCL_NodeGroup::ckNew();

//     Strategy *strategy = new StreamingStrategy(100, 1024, 32*1024, 32*1024*1024);
//     streaming_strat = ComlibRegister(strategy);
//     ComlibAssociateProxy(streaming_strat, VanillaFMMProxy);

    // Get input args
    if (msg->argc < 2) {
        std::cout << "Must specify a file containing charge locations in x,y,z,q\\n format" << std::endl;
        throw std::exception();
    }
    
    std::string filename(msg->argv[1]);
    std::cout << "Hello.  This is the FMM test program.  Reading test data from " << filename << std::endl;
    
    unsigned int max_objects_per_cell = default_max_objects;
    if (msg->argc >= 3) {
        max_objects_per_cell = boost::lexical_cast<unsigned int>(msg->argv[2]);
        std::cout << "Using user-specified max_objects_per_cell=" << max_objects_per_cell << std::endl;
    
    }

    srand((unsigned)time(0));
    std::cout.precision(dbl::digits10);
    
    std::ifstream fin;
    fin.open(filename.c_str());
    
    // check that the file stream is valid
    assert(fin.good());

    // list of charges
    std::vector<Charge> charges;

    Vector max(-1e99,-1e99,-1e99);
    Vector min(+1e99,+1e99,+1e99);
    double charge;
    Vector xyz;
    while(fin >> xyz.x >> xyz.y >> xyz.z >> charge)
    {
        // Now remove the extra stuff on the line you do not want.
        fin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
        
        max.x = (xyz.x > max.x) ? xyz.x : max.x;
        max.y = (xyz.y > max.y) ? xyz.y : max.y;
        max.z = (xyz.z > max.z) ? xyz.z : max.z;
        min.x = (xyz.x < min.x) ? xyz.x : min.x;
        min.y = (xyz.y < min.y) ? xyz.y : min.y;
        min.z = (xyz.z < min.z) ? xyz.z : min.z;

        charges.push_back(Charge(xyz,charge));
    }

    // close file
    fin.close();

    const Vector centre = (max + min) / 2.0;
    universe_edge_length = (max - min).length() * 1.1;
    std::cout << "Num Charges: " << charges.size() << " Centre: " << centre << " Edge length: " << universe_edge_length << std::endl;

#if 0
    int ctr=0;
    for(std::vector<Charge>::iterator it1=charges.begin(), end1=charges.end(); it1 != end1; ++it1)
    {
        double pot = 0.0;
        const Vector& pt1 = *it1;
        for(std::vector<Charge>::iterator it2=charges.begin(), end2=charges.end(); it2 != end2; ++it2)
        {
            if (it1 == it2) { continue; }
            pot += it2->charge / (pt1 - *it2).length();
        }
        pot *= ONE_OVER_4PI;
        std::cout << pt1 << ": " << pot << "\n";
        if (++ctr == 10) { break; }
    }
#endif
    
    VanillaFMMProxy = CProxy_VanillaFMM<18,18,300> ::ckNew(max_objects_per_cell, centre, universe_edge_length, charges.size(), CkSelfCallback(CkIndex_Main::create_workers()));
    VanillaFMMProxy.insert(charges);
    
    // We are done with msg so delete it.
    delete msg;

}

void Main::create_workers()
{
    CkPrintf("Solving...\n");
    VanillaFMMProxy.solve_and_evaluate(KAPPA, CkSelfCallback(CkIndex_Main::quiescenceHandler()));
    
}

void Main::fmm_worker_complete()
{
    CkPrintf("Evaluating...\n");
    VanillaFMMProxy.evaluate(CkSelfCallback(CkIndex_Main::quiescenceHandler()));


}
void Main::quiescenceHandler() {

    std::cout << "Complete!" << std::endl;
    std::cout << "Total time: " << (myclock() - start_time)/1000. << " ms" << std::endl;
    
    CkExit();
}

// Constructor needed for chare object migration (ignore for now)
// NOTE: This constructor does not need to appear in the ".ci" file
Main::Main(CkMigrateMessage* msg) {
    std::cout << "AARGH!  Main is migrating!" << std::endl;
    throw std::exception();
}

#define CK_TEMPLATES_ONLY
#include "vanilla_fmm.def.h"
#undef CK_TEMPLATES_ONLY

#include "main.def.h"
