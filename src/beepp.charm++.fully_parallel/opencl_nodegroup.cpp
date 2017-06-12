/*
* opencl_nodegroup.cpp
*
*  Created on: 6 Sep 2010
*      Author: david
*/

#include "prerequisites.h"
#include "opencl_nodegroup.decl.h"
#include "opencl_nodegroup.h"

#include "parallel_fmm_octree.decl.h"
#include "parallel_fmm_octree.h"
#include "fh_values_nodegroup.decl.h"
#include "fh_values_nodegroup.h"
#include "fmm_globals_nodegroup.decl.h"
#include "fmm_worker.decl.h"
#include "fmm_worker.h"
#include "iteration_handler.decl.h"
#include "iteration_handler.h"

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/thread/thread.hpp>

#ifdef OPENCL
#include "../opencl/opencl_handler.h"
#include "../bem/opencl_bem.h"
#endif

extern /* readonly */ beepp::CProxy_ParallelTree ParallelFMMOctreeProxy;
extern /* readonly */ CProxy_FH_Values_NodeGroup FH_Values_NodeGroupProxy;
extern /* readonly */ CProxy_IterationHandler IterationHandlerProxy;

typedef ParallelFMMOctree<CharmNodePatch, CProxy_FMMWorker> ParallelTree;

OpenCL_NodeGroup::OpenCL_NodeGroup(unsigned int num_patches) : m_ready(false), total_num_patches(num_patches) {
    //std::cout << "Created OpenCL NodeGroup" << std::endl;
#ifdef OPENCL
    lock = CmiCreateLock();
    ocl_handler = new OpenCL_Handler;
    //std::cout << "Created OpenCL NodeGroup: total_num_patches=" <<  total_num_patches << std::endl;
    results.reset(new double[total_num_patches*2]);
    memset(results.get(), 0, sizeof(double)*total_num_patches*2);
#endif
    m_ready = true;
}

OpenCL_NodeGroup::~OpenCL_NodeGroup() {
    assert(m_ready);
#ifdef OPENCL
    CmiDestroyLock(lock);
    delete ocl_handler;
#endif
}

void OpenCL_NodeGroup::run_bem(CkArrayIndexOctreeIndexer idxer, double kappa)
{
    assert(m_ready);
#ifdef OPENCL

    CmiLock(lock);

    // get the corresponding node from the FMM tree
    const ParallelTree& tree = *(ParallelFMMOctreeProxy.ckLocalBranch());
    const ParallelTree::NodeT &node = tree.get_node(idxer);
    typedef ParallelTree::NodeT::ContentList ContentList;
    boost::shared_ptr<ContentList> neighbourhood = tree.get_neighbourhood_contents(node);
    
    FH_Values_NodeGroup& fh_vals_group = *(FH_Values_NodeGroupProxy.ckLocalBranch());
    const double *fvals = fh_vals_group.fvals();
    const double *hvals = fh_vals_group.hvals();
    unsigned int total_num_patches = fh_vals_group.get_num_patches();

    // Get the neighbourlist-- i.e. all patches in the adjacent 26 cubes
    size_t num_neighbours = neighbourhood->size();
    size_t ctr=0;
    size_t chunksize=static_cast<size_t>(ceil(static_cast<double>(BEM_EXPLICIT_CHUNKSIZE) / node.size()));
    chunksize = (num_neighbours > chunksize) ? chunksize : num_neighbours;

    while (ctr < num_neighbours)
    {
        PatchPtrList list_chunk(new PPList);
        list_chunk->reserve(chunksize);

        for (size_t ii=0; ii < chunksize && ctr < num_neighbours; ++ii)
        {
            list_chunk->push_back((*neighbourhood)[ctr++]);
        }

        // PatchPtrList is a boost::shared_ptr< std::vector<BasicNodePatch*> > -- i.e. reference
        // counted std::vector of raw node-patch pointers- which should not go out of scope as
        // they reside in the FMM tree.  (Note could make this safer by making the NodePatch* 
        // reference counted, but that would reduce the efficiency of the FMM octree)
        PatchPtrList contents(new PPList);
        contents->insert(contents->begin(), node.get_contents().begin(), node.get_contents().end());

        // IMPORTANT: these heap allocated resources will be deleted by the OpenCL_Handler when they complete execution
        // The memory pointed to be explicit_integrations must remain available, valid, and not be written to by any other
        // threads until the OpenCL_Handler has finished with it.
        BEM_OnDemand_Resources* res_ptr = new BEM_OnDemand_Resources(contents, list_chunk, kappa, fvals, hvals, &(results[0]), &(results[total_num_patches]), true, NULL);
        assert(res_ptr != NULL);
        ocl_handler->add_work_to_queue(res_ptr);
    }
    
    IterationHandlerProxy.done_evals();
    CmiUnlock(lock);

#endif
}

void OpenCL_NodeGroup::precalc_bem(CkArrayIndexOctreeIndexer idxer, double kappa)
{
#ifdef OPENCL

    CmiLock(lock);

    // get the corresponding node from the FMM tree
    const ParallelTree& tree = *(ParallelFMMOctreeProxy.ckLocalBranch());
    const ParallelTree::NodeT &node = tree.get_node(idxer);
    typedef ParallelTree::NodeT::ContentList ContentList;
    boost::shared_ptr<ContentList> neighbourhood = tree.get_neighbourhood_contents(node);
    
    FH_Values_NodeGroup& fh_vals_group = *(FH_Values_NodeGroupProxy.ckLocalBranch());
    const double *fvals = fh_vals_group.fvals();
    const double *hvals = fh_vals_group.hvals();
    unsigned int total_num_patches = fh_vals_group.get_num_patches();
    
    // Get the neighbourlist-- i.e. all patches in the adjacent 26 cubes
    size_t num_neighbours = neighbourhood->size();
    size_t ctr=0;
    size_t chunksize=static_cast<size_t>(ceil(static_cast<double>(BEM_EXPLICIT_CHUNKSIZE) / node.size()));
    chunksize = (num_neighbours > chunksize) ? chunksize : num_neighbours;

    while (ctr < num_neighbours)
    {
        PatchPtrList list_chunk(new PPList);
        list_chunk->reserve(chunksize);

        for (size_t ii=0; ii < chunksize && ctr < num_neighbours; ++ii)
        {
            list_chunk->push_back((*neighbourhood)[ctr++]);
        }

        size_t size = list_chunk->size() * node.size();
        LintArray results(new LocalIntegrations[size]);
        LintArray_Size arr_siz(results, size);
        explicit_integrations.push_back( arr_siz );
        
        // PatchPtrList is a boost::shared_ptr< std::vector<BasicNodePatch*> > -- i.e. reference
        // counted std::vector of raw node-patch pointers- which should not go out of scope as
        // they reside in the FMM tree.  (Note could make this safer by making the NodePatch* 
        // reference counted, but that would reduce the efficiency of the FMM octree)
        PatchPtrList contents(new PPList);
        contents->insert(contents->begin(), node.get_contents().begin(), node.get_contents().end());
        
        // IMPORTANT: these heap allocated resources will be deleted by the OpenCL_Handler when they complete execution
        // The memory pointed to be explicit_integrations must remain available, valid, and not be written to by any other
        // threads until the OpenCL_Handler has finished with it.
        BEM_Resources* res_ptr = new BEM_Resources(contents, list_chunk, kappa, results.get());
        assert(res_ptr != NULL);
        ocl_handler->add_work_to_queue(res_ptr);
    }
    
    
    CmiUnlock(lock);
#endif
}

void OpenCL_NodeGroup::collate_bem_results(CkCallback cb)
{
    assert(m_ready);
#ifdef OPENCL

    CmiLock(lock);

    // if there's still stuff waiting to complete in the opencl handler, then re-call
    int n = static_cast<int>(ocl_handler->pending());
    if (n > 0) {
        CmiUnlock(lock);
        
        CcdCallOnCondition(CcdPERIODIC, resume_thread, CthSelf());
        CthSuspend();
        
        CkEntryOptions *opts = new CkEntryOptions; 
        opts->setPriority(20);
        thisProxy[CkMyNode()].collate_bem_results(cb, opts);
        return;
    }
    
    FH_Values_NodeGroup& fh_vals_group = *(FH_Values_NodeGroupProxy.ckLocalBranch());

    // Run through explicit integrations (if any were precalculated)
    for (std::vector<LintArray_Size>::const_iterator it=explicit_integrations.begin(), end=explicit_integrations.end(); it != end; ++it)
    {
        const LintArray& array = it->first;
        const size_t& sz = it->second;
        for (size_t ii=0; ii < sz; ++ii)
        {
            array[ii].evaluate_local_contributions(fh_vals_group.fvals(), fh_vals_group.hvals(), results.get(), &(results[total_num_patches]));
        }
    }
    
    fh_vals_group.add_results(results.get());
    memset(results.get(), 0, sizeof(double)*total_num_patches*2);

    contribute(cb);

    CmiUnlock(lock);

#else

    // this is a reduction call - must contribute or we'll hang
    contribute(cb);
    
#endif
}

long OpenCL_NodeGroup::add_to_queue(OpenCL_WorkBlob* ptr, bool track)
{
    assert(m_ready);
    long retval = 0;
#ifdef OPENCL
    CmiLock(lock);
    retval = ocl_handler->add_work_to_queue(ptr, track);
    CmiUnlock(lock);
#endif
    return retval;
}

long OpenCL_NodeGroup::add_to_queue_blocking(OpenCL_WorkBlob* ptr, bool track)
{
    assert(m_ready);
    long retval = 0;
#ifdef OPENCL
    CmiLock(lock);
    retval = ocl_handler->add_work_to_queue_blocking(ptr, track);
    CmiUnlock(lock);
#endif
    return retval;
}

void OpenCL_NodeGroup::wait_for_completion()
{
    assert(m_ready);
#ifdef OPENCL
    // this will block until the OpenCL handler is empty
    CmiLock(lock);
    ocl_handler->wait_until_idle();
    CmiUnlock(lock);
#endif
    return;
}

int OpenCL_NodeGroup::pending() {

    assert(m_ready);
#ifdef OPENCL
    CmiLock(lock);
    int pending = static_cast<int>(ocl_handler->pending());
    CmiUnlock(lock);
    return pending;
#else
    return 0;
#endif

}

bool OpenCL_NodeGroup::item_has_finished(long id, bool remove_if_done) 
{
    assert(m_ready);
    bool retval = true;
#ifdef OPENCL
    CmiLock(lock);
    retval = ocl_handler->item_has_finished(id, remove_if_done);
    CmiUnlock(lock);
#endif
    return retval;

}


#include "opencl_nodegroup.def.h"


