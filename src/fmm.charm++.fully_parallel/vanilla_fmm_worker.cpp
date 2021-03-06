#include "vanilla_fmm_prerequisites.h"

#include "fmm_globals_nodegroup.decl.h"
#include "parallel_fmm_octree.decl.h"
#include "vanilla_fmm_worker.decl.h"
#include "vanilla_fmm_evals.decl.h"
#include "vanilla_fmm_worker.h"
#include "vanilla_fmm_evals.h"

#ifdef OPENCL
#include "../fmm/opencl_fmm.h"
#include "opencl_nodegroup.decl.h"
#include "opencl_nodegroup.h"
extern /* readonly */ CProxy_OpenCL_NodeGroup OpenCL_NodeGroupProxy;
#endif

namespace vanilla_fmm
{
    
template<int NTERMS, int NLAMBS, int NWAVES>
inline void VanillaFMMWorkerT<NTERMS,NLAMBS,NWAVES>::evaluate(Eval_Message* msg)
{
    // should only call this when everything is solved... if not ready to evaluate, recall
    // with a low priority
    if (ready_to_evaluate == false)
    {
        // can suspend thread, if this is a threaded function... it's not though
        //CcdCallOnCondition(CcdPERIODIC_1second, resume_thread, CthSelf());
        //CthSuspend();
        super::thisProxy[super::thisIndex].evaluate(msg);
        return;
    }
    
    //std::cout << OctreeIndexer(super::thisIndex) << " : evaluate called! " << std::endl;

    EvalPt* eval_pts = msg->data;
    size_t num_eval_pts = msg->length;
    for (size_t ii=0; ii < num_eval_pts; ++ii)
    {
        eval_pts[ii].reinit_mutex();
    }
    
    const ParallelTree& tree = *(FMM_Tree_Proxy.ckLocalBranch());
    const NodeT& node = tree.get_node(super::thisIndex);

    // sanity check the input -- if nothing to do, there's nothing to do!
    assert(num_eval_pts > 0);

    std::vector<long> tracking_ids;

    typedef boost::shared_ptr< std::vector< Charge* > > NBListPtr;
    
#ifdef OPENCL

    // OpenCL Max memory constraint -- to do the FMM explicit interactions in OpenCL
    // we need to allocate a chunk of global memory for thread-block level results.
    // For large FMM problems it is quite possible that we will exceed the maximum malloc
    // size (e.g. 255MB on GTX280) so need to break the problem into chunks if we find
    // that we are approaching the limit.
    OpenCL_NodeGroup& ocl_nodegroup = *(OpenCL_NodeGroupProxy.ckLocalBranch());
    const OpenCL_Handler& ocl_handler = ocl_nodegroup.get_ocl_handler();
    size_t opencl_max_malloc = ocl_handler.get_max_malloc_size();

    NBListPtr neighbourhood = tree.get_neighbourhood_contents(node);

    // break problem into smaller chunks if necessary
    unsigned int num_x_chunks = FMM_Resources<EvalPt>::calc_x_chunksize(num_eval_pts, neighbourhood->size(), opencl_max_malloc);
    unsigned int num_y_chunks = FMM_Resources<EvalPt>::calc_y_chunksize(num_eval_pts, neighbourhood->size(), opencl_max_malloc);

    // if number of chunks is 1, then no need to split the problem up
    if (num_x_chunks == 1 && num_y_chunks == 1)
    {
        // IMPORTANT: this is heap allocated and then passed to the OpenCL handler which process the
        // work in a separate thread-- it is the job of the OpenCL handler to correctly delete this!
        std::vector<EvalPt*> eval_pts_chunk;
        eval_pts_chunk.reserve(num_eval_pts);
        for (size_t ii=0; ii < num_eval_pts; ++ii)
        {
            eval_pts_chunk.push_back(eval_pts + ii);
        }
        FMM_Resources<EvalPt>* fmm_res_ptr = new FMM_Resources<EvalPt>(eval_pts_chunk, *neighbourhood, beta);
        long opencl_tracking_id = ocl_nodegroup.add_to_queue(fmm_res_ptr, true);
        tracking_ids.push_back(opencl_tracking_id);
    }
    else
    {
        // Annoyingly complicated: break the neighbourhood into equal sized lumps, and
        // create an FMM work unit for each
        size_t x_chunk_size = static_cast<size_t>(ceil(static_cast<double>(num_eval_pts) / num_x_chunks));
        size_t y_chunk_size = static_cast<size_t>(ceil(static_cast<double>(neighbourhood->size()) / num_y_chunks));

        size_t eval_ctr=0;
        while(eval_ctr != num_eval_pts)
        {
            std::vector<EvalPt*> eval_pts_chunk;
            eval_pts_chunk.reserve(x_chunk_size);
            for (size_t x_chunk_ctr=0 ; eval_ctr != num_eval_pts && x_chunk_ctr < x_chunk_size; ++eval_ctr)
            {
                eval_pts_chunk.push_back(eval_pts + eval_ctr);
                ++x_chunk_ctr;
            }

            typename std::vector< Charge* >::const_iterator neigh_it=neighbourhood->begin(), neigh_end=neighbourhood->end();
            while(neigh_it != neigh_end)
            {
                std::vector< Charge* > neighbourhood_chunk;
                neighbourhood_chunk.reserve(y_chunk_size);
                for (size_t y_chunk_ctr=0 ; neigh_it != neigh_end && y_chunk_ctr < y_chunk_size; ++neigh_it)
                {
                    neighbourhood_chunk.push_back(*neigh_it);
                    ++y_chunk_ctr;
                }
                FMM_Resources<EvalPt>* fmm_res_ptr = new FMM_Resources<EvalPt>(eval_pts_chunk, neighbourhood_chunk, beta);
                long opencl_tracking_id = ocl_nodegroup.add_to_queue(fmm_res_ptr, true);
                tracking_ids.push_back(opencl_tracking_id);

            }
        }
    }
#else
    NBListPtr neighbourhood = tree.get_neighbourhood_contents(node);
    for (size_t ii=0; ii < num_eval_pts; ++ii)
    {
        for (std::vector< Charge* >::const_iterator neigh_it=neighbourhood->begin(), neigh_end=neighbourhood->end(); neigh_it != neigh_end; ++neigh_it)
        {
            EvalPt::add_explicit_contrib(eval_pts[ii], **neigh_it, beta);
        }
    }
#endif

    unsigned short level = node.get_level();
    const double scale = FMM_Globals<NTERMS>::get_scale(beta, tree.get_edge_length(), level);
    const Vector& node_centre = node.get_centre();
    const MultipoleHolder& local_expansions = *lpole_ptr;

    if (level >= 2) {

        // loop over eval points
        for (size_t ii=0; ii < num_eval_pts; ++ii)
        {
            double pot=0;
            Vector field(0,0,0);
            fmm::evaluate_local_expansion_at_xyz(beta,
                                                 scale,
                                                 eval_pts[ii],
                                                 node_centre,
                                                 local_expansions);
        }
    }
/*
#ifdef OPENCL
    while (tracking_ids.size())
    {
        // iterate over the list of tracking ids (the work might've been split into multiple chunks) and 
        // remove items which have completed, if they're all completed then hooray job done
        std::vector< std::vector<long>::iterator > erase_list;
        for (std::vector<long>::iterator track_it=tracking_ids.begin(), track_end=tracking_ids.end(); track_it != track_end; ++track_it)
        {
            // this will remove the tracking_id from the opencl handler if complete- so if it returns true, don't ask about the same id again!
            if (ocl_nodegroup.item_has_finished(*track_it, true)) { // 'true' here enforces removal of complete items
                //std::cout << "Tracking item: " << *track_it << " completed." << std::endl;
                erase_list.push_back(track_it); // don't ask multiple times for same item 
            }
        }
        
        // erase the erase-list (i wish the stl containers could erase with a reverse_iterator...)
        // this is well cumbersome :-(
        for (std::vector< std::vector<long>::iterator >::reverse_iterator erase_it=erase_list.rbegin(), erase_end=erase_list.rend(); erase_it!=erase_end; ++erase_it) {
            tracking_ids.erase(*erase_it);
        }
        erase_list.clear();
        
        if (tracking_ids.size())
        {
            // suspend this thread
            CcdCallOnCondition(CcdPERIODIC_100ms, resume_thread, CthSelf());
            CthSuspend();
        }
    }
    
    if (tracking_idx.size() == 0)
    {
        msg->cb.send(msg);
        return;
    }

#endif*/

#ifdef OPENCL
    if (tracking_ids.size() > 0)
    {
         work_pending.push_back(Msg_Tracking_Pair(msg, tracking_ids));
         
         // only actually need to call the asynchronous check once and
         // then it'll re-call itself until the work_pending list is empty, 
         // so make the call only if the thing we just added was the first
         // one in the list.  Note that since charm++ only allows one entry
         // method per chare to run there's no chance of dodgy races and thready
         // stuff happening with this bit, so this should work fine.
         if (work_pending.size() == 1)
         {
            super::thisProxy[super::thisIndex].asynchronous_check();
         }
    }
    else
    {
        //std::cout << OctreeIndexer(super::thisIndex) << " : Completed work, sending message..." << std::endl;
        msg->cb.send(msg);
    }
    
#else

    msg->cb.send(msg);

#endif

}

template<int NTERMS, int NLAMBS, int NWAVES>
void VanillaFMMWorkerT<NTERMS,NLAMBS,NWAVES>::asynchronous_check(int dummy)
{
    
#ifdef OPENCL
    while(true)
    {

        // nothing to do if the work_pending list is empty
        if (work_pending.size() == 0) { return; }
        
        // get the local OpenCL handler
        OpenCL_NodeGroup& ocl_nodegroup = *(OpenCL_NodeGroupProxy.ckLocalBranch());
        
        // this checks the list of OpenCL tracking ids (see completed_tracked list in OpenCL handler) to figure out
        // whether the OpenCL bits of evaluation are completed; if so then the corresponding message can be sent
        // through the embedded callback.  Phew- is all a bit complicated...
        std::vector<std::vector<Msg_Tracking_Pair>::iterator> work_pending_erase_list;
        for (std::vector<Msg_Tracking_Pair>::iterator it=work_pending.begin(), end=work_pending.end(); it != end; ++it)
        {
            Eval_Message* msg = it->first;
            std::vector<long>& tracking_ids = it->second;
            
            // iterate over the list of tracking ids (the work might've been split into multiple chunks) and 
            // remove items which have complete, if they're all completed then hooray job done
            std::vector< std::vector<long>::iterator > erase_list;
            for (std::vector<long>::iterator track_it=tracking_ids.begin(), track_end=tracking_ids.end(); track_it != track_end; ++track_it)
            {
                // this will remove the tracking_id from the opencl handler if complete- so if it returns true, don't ask about the same id again!
                if (ocl_nodegroup.item_has_finished(*track_it, true)) { // 'true' here enforces removal of complete items
                    //std::cout << "Tracking item: " << *track_it << " completed." << std::endl;
                    erase_list.push_back(track_it); // don't ask multiple times for same item 
                }
            }
            
            // erase the erase-list (i wish the stl containers could erase with a reverse_iterator...)
            // this is well cumbersome :-(
            for (std::vector< std::vector<long>::iterator >::iterator erase_it=erase_list.begin(), erase_end=erase_list.end(); erase_it!=erase_end; ++erase_it) {
                tracking_ids.erase(*erase_it);
            }

            if (tracking_ids.size() == 0) // if not waiting for any more OpenCL chunks, then all are complete
            {
                // send the message, it is complete
                //std::cout << OctreeIndexer(super::thisIndex) << " : [async] Completed work, sending message..." << std::endl;
                msg->cb.send(msg);
                
                // remove the Msg_Tracking_Pair
                work_pending_erase_list.push_back(it);
            }

        }
        
        // erase the erase-list -- this time for the work_pending list (it's still well cumbersome)
        for (std::vector<std::vector<Msg_Tracking_Pair>::iterator>::iterator erase_it=work_pending_erase_list.begin(), erase_end=work_pending_erase_list.end(); erase_it!=erase_end; ++erase_it) {
            work_pending.erase(*erase_it);
        }
        
        // if there's still work pending, then need to re-trigger this function in the future
        int n = work_pending.size();
        if (n > 0)
        {
            CkEntryOptions* opts(new CkEntryOptions);
            opts->setPriority(20);
            super::thisProxy[super::thisIndex].asynchronous_check(0, opts);
        }
        else
        {
            break;
        }
    }

#endif

}

} // end namespace

#define CK_TEMPLATES_ONLY
#include "parallel_fmm_octree.def.h"
#undef CK_TEMPLATES_ONLY

#include "vanilla_fmm_worker.def.h"
