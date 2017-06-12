// things that appear in the .decl files and therefore need to be included first
#include "vanilla_fmm_prerequisites.h"

// Other Chare types
#include "main.decl.h"
#include "fmm_globals_nodegroup.decl.h"
#include "parallel_fmm_octree.decl.h"
#include "parallel_fmm_octree.h"
#include "vanilla_fmm_worker.decl.h"
#include "vanilla_fmm_worker.h"
#include "vanilla_fmm_evals.decl.h"
#include "vanilla_fmm_evals.h"
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

typedef ParallelFMMOctree<Charge, vanilla_fmm::CProxy_VanillaFMMWorker> ParallelTree;

/* readonly */ CProxy_Main MainProxy;
/* readonly */ CProxy_ParallelFMMOctree<Charge, vanilla_fmm::CProxy_VanillaFMMWorker> ParallelFMMOctreeProxy;
/* readonly */ CProxy_OpenCL_NodeGroup OpenCL_NodeGroupProxy;
/* readonly */ vanilla_fmm::CProxy_FMM_Globals_NodeGroup FMM_Globals_Proxy;
/* readonly */ vanilla_fmm::CProxy_VanillaFMMWorker VanillaFMMWorkerProxy;
/* readonly */ vanilla_fmm::CProxy_Vanilla_FMM_Evals Vanilla_FMM_Evals_Proxy;
/* readonly */ ComlibInstanceHandle streaming_strat;

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
Main::Main(CkArgMsg* msg) : pending(0)
{

    CkPrintf("Running Parallel FMM using %d processors.\n", CkNumPes());

    // Set the mainProxy readonly to point to a
    // proxy for the Main chare object (this
    // chare object).
    MainProxy = thisProxy;
    FMM_Globals_Proxy = CProxy_FMM_Globals_NodeGroup::ckNew();
    VanillaFMMWorkerProxy = CProxy_VanillaFMMWorker::ckNew();
    OpenCL_NodeGroupProxy = CProxy_OpenCL_NodeGroup::ckNew();

    Strategy *strategy = new StreamingStrategy(100, 1024, 32*1024, 32*1024*1024);
    streaming_strat = ComlibRegister(strategy);
    ComlibAssociateProxy(streaming_strat, ParallelFMMOctreeProxy);
    //ComlibAssociateProxy(streaming_strat, VanillaFMMWorkerProxy);

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
        eval_pts.push_back( xyz ); // this is stored in the Main class- need it later to create the Vanilla_FMM_Evals chare
    }

    // close file
    fin.close();

    const Vector centre = (max + min) / 2.0;
    universe_edge_length = (max - min).length() * 1.1;
    const unsigned int num_eval_pts = charges.size();
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
            double r = (pt1-*it2).length();
            pot += it2->charge *exp(-KAPPA*r) / r;
        }
        //pot *= ONE_OVER_4PI;
        std::cout << pt1 << ": " << pot << "\n";
        if (++ctr == 10) { break; }
    }
#endif
    
    start_time = myclock();
    ParallelFMMOctreeProxy = CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker>::ckNew(max_objects_per_cell, centre, universe_edge_length, charges.size(), CkSelfCallback(CkIndex_Main::create_workers()));
    ParallelFMMOctreeProxy.insert(charges);
    
    // We are done with msg so delete it.
    delete msg;

}

void Main::create_workers()
{
    
    CkPrintf("Creating evaluators...\n");
    Vanilla_FMM_Evals_Proxy = CProxy_Vanilla_FMM_Evals::ckNew(eval_pts, ParallelFMMOctreeProxy, VanillaFMMWorkerProxy);

    CkPrintf("Creating workers...\n");
    unsigned int num_workers=0, total_nodes=0;
    
    // Create worker chares
    ParallelTree *parallel_octree = ParallelFMMOctreeProxy.ckLocalBranch();
    assert(parallel_octree != NULL);
    for (unsigned short level=parallel_octree->get_top_level(); level <= parallel_octree->get_bottom_level(); ++level)
    {
        for (ParallelTree::NodeList::const_iterator it=parallel_octree->get_node_list(level).begin(), end=parallel_octree->get_node_list(level).end();
                it != end;
                ++it)
        {
            const ParallelTree::NodeT& node = *(it->second);
            assert(node.empty() == false);
            VanillaFMMWorkerProxy[node.get_idx()].insert(FMM_Globals_Proxy, ParallelFMMOctreeProxy, universe_edge_length);
            if (node.isLeaf()) {
                
                ++num_workers;
            }
            ++total_nodes;
        }
    }
    VanillaFMMWorkerProxy.doneInserting();
    std::cout << "info: " << total_nodes << " chares, of which " << num_workers << " are leaves" << std::endl;
    pending = num_workers; // number of callback events to expect when solving FMM
    
    CkPrintf("Solving...\n");
    VanillaFMMWorkerProxy.solve(KAPPA, CkSelfCallback(CkIndex_Main::fmm_worker_complete()));
    
    // Can allow load balancing now
    //TurnManualLBOff();
    //LBDatabaseObj()->StartLB();
    
}

void Main::fmm_worker_complete()
{
    if(--pending == 0)
    {
        std::cout << "Evaluating..." << std::endl;        
        
        // Create the FMM Evaluation proxy- this controls the actual evaluations: need FMM to complete before can
        // do FMM evals (obviously), but explicit work can be done sooner- can backfill some of the idle CPU cycles.
        Vanilla_FMM_Evals_Proxy.evaluate(CkSelfCallback(CkIndex_Main::completed(NULL)));
        
    }

}

void Main::completed(Eval_Message* msg)
{
    std::cout << "Finished!" << std::endl;
    
    size_t num_pts = msg->length;
    EvalPt* eval_pts = msg->data;
    for (size_t ii=0; ii < num_pts && ii < 10; ++ii)
    {
        const EvalPt& ep = eval_pts[ii];
        std::cout << ep.pt() << ": " << ep.get_potential() << std::endl;
    }
    delete msg;
 
    std::cout << "Total time: " << (myclock() - start_time)/1000. << " ms" << std::endl;
    
    CkExit();
    //CkStartQD(CkCallback(CkIndex_Main::quiescenceHandler(), thisProxy));
}

void Main::quiescenceHandler() {

    std::cout << "Quiescence!" << std::endl;
    eval_pts.clear();
    
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
#include "parallel_fmm_octree.def.h"
//#include "fmm_globals_nodegroup.def.h"
#undef CK_TEMPLATES_ONLY

#include "main.def.h"
