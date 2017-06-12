// things that appear in the .decl files and therefore need to be included first
// basically a bunch of C++ headers defining the various BEM/FMM objects
#include "prerequisites.h"

// Other Chare types
#include "fh_values_nodegroup.decl.h"
#include "mesh_working_copy.decl.h"
#include "mesh_library.decl.h"
#include "mesh_library.h"
#include "rhs_handler.decl.h"
#include "rhs_handler.h"
#include "gmres.decl.h"
#include "gmres.h"
#include "fh_values_nodegroup.decl.h"
#include "iteration_handler.decl.h"
#include "../bem/local_integrations.h"
#include "bem_fmm.decl.h"
#include "bem_fmm.h"

#include "main.decl.h"
#include "main.h"

// 'Normal' C++ headers
#include <string>
#include <boost/program_options.hpp>

namespace po = boost::program_options;
/* readonly */ CProxy_BEM_FMM<NT,NL,NW,CharmNodePatch> BEM_FMM_Proxy;

/* readonly */ CProxy_Main MainProxy;
/* readonly */ CProxy_MeshWorkingCopy MeshWorkingCopyProxy;
/* readonly */ CProxy_MeshLibrary MeshLibraryProxy;
/* readonly */ CProxy_RHS_Handler RHS_HandlerProxy;
/* readonly */ CProxy_GMRES GMRESProxy;
/* readonly */ CProxy_IterationHandler IterationHandlerProxy;
/* readonly */ CProxy_FH_Values_NodeGroup FH_Values_NodeGroupProxy;
/* readonly */ CProxy_OpenCL_NodeGroup OpenCL_NodeGroupProxy;

// Entry point of Charm++ application
Main::Main(CkArgMsg* msg) : default_quad_points_per_triangle(DEFAULT_QUADS), default_qual_points_per_triangle(DEFAULT_QUALS)
{

    start_clock = myclock();
    
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
    po::store(po::command_line_parser(msg->argc, msg->argv).
            options(opts).positional(p).run(), vm);
    po::notify(vm);
    
    // Read the config file which describes the simulation set-up
    //std::string config_filename(msg->argv[1]);

    // We are done with msg so delete it.
    delete msg;
    
    CkPrintf("Running BEM/FMM using %d processors.\n", CkNumPes());
    //CkPrintf("Sizeof CProxy_MeshLibrary: %d bytes.\n", sizeof(CProxy_MeshLibrary));
    //CkPrintf("Sizeof CharmNodePatch: %d bytes.\n", sizeof(CharmNodePatch));

    // set global readonly vars (proxies to parallel objects)
    MainProxy = thishandle;
    GMRESProxy = CProxy_GMRES::ckNew();
    
    // Create ConfigFile from the xml config file
    config.reset(new ConfigFile(config_filename));
    output_filename = config->output_file;

    Dsolvent = config->solvent_dielectric;
    kappa = (cmdline_kappa == -1) ? config->solvent_kappa : cmdline_kappa;
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

        std::cout <<"default Quads per triangle: " << default_quad_points_per_triangle << std::endl;
        std::cout <<"default Quals per triangle: " << default_qual_points_per_triangle << std::endl;

        std::cout <<"xml Quads per triangle: " << it->quad_points << std::endl;
        std::cout <<"xml Quals per triangle: " << it->qual_points << std::endl;
        
        // respect the settings for quad/qual points in the xml configuration...
        unsigned int quad_points_per_triangle = (it->quad_points == -1) ? default_quad_points_per_triangle : static_cast<unsigned int>(it->quad_points);
        unsigned int qual_points_per_triangle = (it->qual_points == -1) ? default_qual_points_per_triangle : static_cast<unsigned int>(it->qual_points);
        
        std::cout <<"Quads per triangle: " << quad_points_per_triangle << std::endl;
        std::cout <<"Quals per triangle: " << qual_points_per_triangle << std::endl;
        
        // ... but global settings on the command line trump everything
        if (cmdline_quad_points != -1) quad_points_per_triangle = static_cast<unsigned int>(cmdline_quad_points);
        if (cmdline_qual_points != -1) qual_points_per_triangle = static_cast<unsigned int>(cmdline_qual_points);

        //std::cout << "Init'ing working mesh id: " << mesh_instance_id << " (" << mesh_lib_id << " " << rotation << " " << xyz_offset << ")" << std::endl;
        MeshWorkingCopyProxy[mesh_instance_id].init(mesh_lib_id, rotation, xyz_offset, total_num_patches, it->dielectric, Dsolvent, quad_points_per_triangle, qual_points_per_triangle);
        total_num_patches += mesh_lib.get_mesh(mesh_lib_id).get_num_node_patches();
        total_num_charges += mesh_lib.get_mesh(mesh_lib_id).get_num_charges();
        
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
    Vector diff = (max - min); 
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
    if (cmdline_nbsize != -1) bem_size = static_cast<unsigned int>(cmdline_nbsize);
    CkPrintf("BEM neighbourhood size: %d\n", bem_size);
    
    // These two octree-related things need the universe edge length and centre
    // ParallelFMMOctreeProxy is the BEM/FMM hybrid octree containing CharmNodePatch objects
    // RHS_HandlerProxy wraps another FMM octree, holding Charge objects, which is used for
    // calculating the RHS vector (using the 'vanilla' fmm rather than the hybrid bem/fmm)
    BEM_FMM_Proxy = CProxy_BEM_FMM<NT,NL,NW,CharmNodePatch>::ckNew(bem_size, universe_centre, universe_edge_length, total_num_patches, CkSelfCallback(CkIndex_Main::get_rhs()));
    RHS_HandlerProxy = CProxy_RHS_Handler::ckNew(MAX_FMM_SIZE, universe_centre, universe_edge_length, total_num_charges);

    FH_Values_NodeGroupProxy = CProxy_FH_Values_NodeGroup::ckNew(total_num_patches);
    OpenCL_NodeGroupProxy = CProxy_OpenCL_NodeGroup::ckNew();
    IterationHandlerProxy = CProxy_IterationHandler::ckNew(universe_edge_length, kappa);


}

// this is a threaded method, which will suspend until we get the RHS data
void Main::get_rhs()
{
    // DEBUG
    CkPrintf("Done creating molecules in simulation space.\n");
    
    CkPrintf("Starting rhs calcs.\n");
    RHS_Message* rhs_msg;
    RHS_HandlerProxy.get_rhs(Dsolvent, total_num_patches, CkCallbackResumeThread((void*&)rhs_msg));
    
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
#include "bem_fmm.def.h"
#include "vanilla_fmm.def.h"
#undef CK_TEMPLATES_ONLY

#include "main.def.h"
