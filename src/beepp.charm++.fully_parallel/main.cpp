// things that appear in the .decl files and therefore need to be included first
// basically a bunch of C++ headers defining the various BEM/FMM objects
#include "prerequisites.h"

// Other Chare types
#include "fh_values_nodegroup.decl.h"
#include "parallel_fmm_octree.decl.h"
#include "parallel_fmm_octree.h"
#include "vanilla_fmm_worker.decl.h"
//#include "vanilla_fmm_worker.h"
#include "vanilla_fmm_evals.decl.h"
#include "mesh_working_copy.decl.h"
#include "mesh_library.decl.h"
#include "mesh_library.h"
#include "fmm_globals_nodegroup.decl.h"
#include "fmm_globals_nodegroup.h"
#include "rhs_handler.decl.h"
#include "rhs_handler.h"
#include "gmres.decl.h"
#include "gmres.h"
#include "explicit_worker.decl.h"
#include "fh_values_nodegroup.decl.h"
#include "iteration_handler.decl.h"
#include "fmm_worker.decl.h"
#include "fmm_worker.h"
#include "../bem/local_integrations.h"

#include "main.decl.h"
#include "main.h"

// 'Normal' C++ headers
#include <string>

using beepp::CProxy_FMM_Globals_NodeGroup;
using beepp::CProxy_FMMWorker;

typedef ParallelFMMOctree<CharmNodePatch, CProxy_FMMWorker> ParallelTree;

/* readonly */ CProxy_Main MainProxy;
/* readonly */ CProxy_MeshWorkingCopy MeshWorkingCopyProxy;
/* readonly */ beepp::CProxy_FMM_Globals_NodeGroup FMM_Globals_Proxy;
/* readonly */ CProxy_MeshLibrary MeshLibraryProxy;
/* readonly */ beepp::CProxy_ParallelTree ParallelFMMOctreeProxy;
/* readonly */ CProxy_FMMWorker FMMWorkerProxy;
/* readonly */ CProxy_ExplicitWorker ExplicitWorkerProxy;
/* readonly */ CProxy_RHS_Handler RHS_HandlerProxy;
/* readonly */ CProxy_GMRES GMRESProxy;
/* readonly */ CProxy_IterationHandler IterationHandlerProxy;
/* readonly */ CProxy_FH_Values_NodeGroup FH_Values_NodeGroupProxy;
/* readonly */ CProxy_OpenCL_NodeGroup OpenCL_NodeGroupProxy;
/* readonly */ ComlibInstanceHandle streaming_strat;

// Entry point of Charm++ application
Main::Main(CkArgMsg* msg)
{
    start_clock = myclock();
    
    CkPrintf("Running BEM/FMM using %d processors.\n", CkNumPes());
    //CkPrintf("Sizeof CProxy_MeshLibrary: %d bytes.\n", sizeof(CProxy_MeshLibrary));
    //CkPrintf("Sizeof CharmNodePatch: %d bytes.\n", sizeof(CharmNodePatch));

    // set global readonly vars (proxies to parallel objects)
    MainProxy = thishandle;
    FMM_Globals_Proxy = beepp::CProxy_FMM_Globals_NodeGroup::ckNew();
    FMMWorkerProxy = CProxy_FMMWorker::ckNew();
    ExplicitWorkerProxy = CProxy_ExplicitWorker::ckNew();
    GMRESProxy = CProxy_GMRES::ckNew();

    Strategy *strategy = new StreamingStrategy(100, 1024, 32*1024, 32*1024*1024);
    streaming_strat = ComlibRegister(strategy);
    ComlibAssociateProxy(streaming_strat, ParallelFMMOctreeProxy);
    
    // Read the config file which describes the simulation set-up
    std::string config_filename(msg->argv[1]);

    // We are done with msg so delete it.
    delete msg;

    // Create ConfigFile from the xml config file
    config.reset(new ConfigFile(config_filename));
    output_filename = config->output_file;

    Dsolvent = config->solvent_dielectric;
    kappa = config->solvent_kappa;
    CkPrintf("Kappa: %f\n", kappa);

    std::vector<std::string> mesh_filenames;
    for (ConfigFile::MeshLibrary::const_iterator it=config->mesh_library.begin(), end=config->mesh_library.end(); it != end; ++it)
    {
        mesh_filenames.push_back(it->mesh_filename);
    }
    MeshLibraryProxy = CProxy_MeshLibrary::ckNew(mesh_filenames);

    const MeshLibrary& mesh_lib = *(MeshLibraryProxy.ckLocalBranch());

    size_t num_working_meshes = config->layout.size();
    MeshWorkingCopyProxy = CProxy_MeshWorkingCopy::ckNew(num_working_meshes);

    // bounding box
    Vector max(-1e99,-1e99,-1e99); // this will be the 'top right' corner
    Vector min(1e99,1e99,1e99); // and this will be 'bottom left'
    
    // loop over the mesh instances in the config layout
    total_num_patches=0;
    for (ConfigFile::Layout::const_iterator it=config->layout.begin(), end=config->layout.end(); it != end; ++it)
    {
        // extract the mesh id
        unsigned int mesh_lib_id = it->mesh_id;
        const Vector& xyz_offset = it->offset;
        const Quaternion& rotation = it->rotation;
        int mesh_instance_id = it->instance_id;

        //std::cout << "Init'ing working mesh id: " << mesh_instance_id << " (" << mesh_lib_id << " " << rotation << " " << offset << ")" << std::endl;
        MeshWorkingCopyProxy[mesh_instance_id].init(mesh_lib_id, rotation, xyz_offset, total_num_patches, it->dielectric, Dsolvent);
        total_num_patches += mesh_lib.get_mesh(mesh_lib_id).get_num_node_patches();
        
        // calculate bounding box for all meshes
        // (need to know the size of the universe in order
        // to init the FMM octrees).
        Vector v;
        double radius = mesh_lib.get_mesh(mesh_lib_id).get_radius();
        v = xyz_offset + Vector(radius, radius, radius);
        max.x = v.x > max.x ? v.x : max.x;
        max.y = v.y > max.y ? v.y : max.y;
        max.z = v.z > max.z ? v.z : max.z;

        min.x = v.x < min.x ? v.x : min.x;
        min.y = v.y < min.y ? v.y : min.y;
        min.z = v.z < min.z ? v.z : min.z;

        v = xyz_offset - Vector(radius, radius, radius);
        max.x = v.x > max.x ? v.x : max.x;
        max.y = v.y > max.y ? v.y : max.y;
        max.z = v.z > max.z ? v.z : max.z;

        min.x = v.x < min.x ? v.x : min.x;
        min.y = v.y < min.y ? v.y : min.y;
        min.z = v.z < min.z ? v.z : min.z;
    }
    
    // figure out the maximum edge length in x/y/z dimension
    Vector diff = (max - min) * 1.02; // add 2% for safety (should not be necessary as radii are maximum boundaries)
    universe_edge_length = diff.y > diff.x ? diff.y : diff.x;
    universe_edge_length = diff.z > universe_edge_length ? diff.z : universe_edge_length;

    assert(universe_edge_length > 0.0);

    // centre of the cube is the mid point of the two extremities
    Vector universe_centre = (max + min) / 2.0;
    
    std::ostringstream buf;
    buf << "Info: universe centre & edge length is: " << universe_centre << " " << universe_edge_length << std::endl;
    buf << "Info: total number of patches is: " << total_num_patches << std::endl;
    CkPrintf(buf.str().c_str());

    // can override the default BEM average neighbourhood from the config filename
    // (add explicit_bem_interactions=XXX parameter to the root xml node)
    unsigned int bem_size = (config->target_bem_explicits != 0) ? config->target_bem_explicits : DEFAULT_BEM_NEIGHBOURHOOD_SIZE;
    CkPrintf("BEM neighbourhood size: %d\n", bem_size);
    
    // These two octree-related things need the universe edge length and centre
    // ParallelFMMOctreeProxy is the BEM/FMM hybrid octree containing CharmNodePatch objects
    // RHS_HandlerProxy wraps another FMM octree, holding Charge objects, which is used for
    // calculating the RHS vector (using the 'vanilla' fmm rather than the hybrid bem/fmm)
    ParallelFMMOctreeProxy = beepp::CProxy_ParallelTree::ckNew(bem_size, universe_centre, universe_edge_length, total_num_patches, CkSelfCallback(CkIndex_Main::get_rhs()));
    RHS_HandlerProxy = CProxy_RHS_Handler::ckNew(MAX_FMM_SIZE, universe_centre, universe_edge_length);

    FH_Values_NodeGroupProxy = CProxy_FH_Values_NodeGroup::ckNew(total_num_patches);
    OpenCL_NodeGroupProxy = CProxy_OpenCL_NodeGroup::ckNew(total_num_patches);
    IterationHandlerProxy = CProxy_IterationHandler::ckNew(universe_edge_length, kappa);


}

// this is a threaded method, which will suspend until we get the RHS data
void Main::get_rhs()
{
    static bool done=false;

    if (done) { return; }
    done = true;

    // DEBUG
    CkPrintf("Done creating molecules in simulation space.\n");

    // create the worker chares
    MainProxy.create_workers(CkCallbackResumeThread());

    // Can allow load balancing now
    LBDatabase * myLbdb = LBDatabase::Object();
    if (myLbdb) {
        myLbdb->TurnManualLBOff();
        myLbdb->StartLB();
    }
    else 
    {
        CkPrintf("Couldn't enable the load balancer.  How annoying.\n");
        //LBDatabase::manualOn = 0;
        //LBDatabase::ckLocalBranch()->StartLB();
    }
    
    CkPrintf("Starting rhs calcs.\n");
    RHS_Message* rhs_msg;
    CkCallbackResumeThread *cb_ptr = new CkCallbackResumeThread((void*&)rhs_msg);
    RHS_HandlerProxy.get_rhs(Dsolvent, total_num_patches, *cb_ptr);
    
    // do the precalcs
    CkPrintf("Carrying out BEM explicit integral pre-calcs in background\n");
    CkEntryOptions *opts = new CkEntryOptions; 
    opts->setPriority(20);
    ExplicitWorkerProxy.precalc(kappa, total_num_patches, opts);
    
    // this'll freeze the thread
    delete cb_ptr;
    
    // send the gmres
    FH_Values rhs_marshallable(rhs_msg->length, rhs_msg->data);
    delete rhs_msg;

    // this is a synchronous call (will suspend this thread)
    GMRESProxy.solve(output_filename, kappa, Dsolvent, total_num_patches, rhs_marshallable);

    // end handler
    quiescenceHandler();

}

void Main::create_workers(CkCallback cb)
{
    CkPrintf("Starting create_workers\n");

    unsigned int num_workers=0;
    // Create worker chares
    ParallelTree *parallel_octree = ParallelFMMOctreeProxy.ckLocalBranch();
    
    unsigned int total_local_ints=0;
    total_local_ints = parallel_octree->calc_neighbourhood_interacts();
    double mem = sizeof(LocalIntegrations)*total_local_ints/(1024.*1024.);
    CkPrintf("Total explicit integrals: %d (%f MB)\n", total_local_ints, mem);
#ifdef CACHE_GPU_RESULTS
    CkPrintf("Explicit integrations *will* be precalculated\n");
#endif
    
    assert(parallel_octree != NULL);
    for (unsigned short level=parallel_octree->get_top_level(); level <= parallel_octree->get_bottom_level(); ++level)
    {
        for (ParallelTree::NodeList::const_iterator it=parallel_octree->get_node_list(level).begin(), end=parallel_octree->get_node_list(level).end();
                it != end;
                ++it)
        {
            const ParallelTree::NodeT& node = *(it->second);
            assert(node.empty() == false);
            FMMWorkerProxy[node.get_idx()].insert(FMM_Globals_Proxy, 
                                                  ParallelFMMOctreeProxy, 
                                                  FH_Values_NodeGroupProxy,
                                                  universe_edge_length);
            if (node.isLeaf())
            {
                ExplicitWorkerProxy[node.get_idx()].insert();
                num_workers += 2;
            }
        }
    }
    FMMWorkerProxy.doneInserting();
    ExplicitWorkerProxy.doneInserting();

    std::ostringstream buf;
    //buf << "hybrid BEM/FMM info: " << num_workers << " worker chares" << std::endl;
    CkPrintf(buf.str().c_str());
    
    IterationHandlerProxy.set_num_workers(num_workers);

    cb.send();

}

void Main::quiescenceHandler() {

    //std::cout << "Quiescence!" << std::endl;
    std::ostringstream buf;
    buf << "Total (wallclock) time (ms): " << (myclock() - start_clock) / 1000. << std::endl;
    CkPrintf(buf.str().c_str());

    CkExit();
}

// Constructor needed for chare object migration (ignore for now)
// NOTE: This constructor does not need to appear in the ".ci" file
Main::Main(CkMigrateMessage* msg) {
    std::cerr << "AARGH!  Main is migrating!" << std::endl;
}

#define CK_TEMPLATES_ONLY
#include "parallel_fmm_octree.def.h"
#undef CK_TEMPLATES_ONLY

#include "main.def.h"
