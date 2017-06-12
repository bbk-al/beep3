#ifndef __VANILLA_FMM_EVALS_H
#define __VANILLA_FMM_EVALS_H

#include "vanilla_fmm_prerequisites.h"
#include "charm++.h"
#include <boost/shared_array.hpp>
#include "../fmm/eval_pt.h"

using fmm::EvalPt;

namespace vanilla_fmm
{

class Eval_Message : public CMessage_Eval_Message {

public:

    int length;
    int length_cb;
    EvalPt* data;
    CkCallback cb;
    OctreeIndexer target;

};

class Vanilla_FMM_Evals : public CBase_Vanilla_FMM_Evals
{
    
    typedef ParallelFMMOctree<Charge, CProxy_VanillaFMMWorker> ParallelTree;
    
public:

    /// Constructors ///
    Vanilla_FMM_Evals(std::vector<Vector> eval_pts_incoming,
                      CProxy_ParallelTree ParallelFMMOctreeProxy, 
                      CProxy_VanillaFMMWorker vanilla_worker_proxy)
    {
        
        VanillaFMMWorkerProxy = vanilla_worker_proxy;
        
        pending = 0;
        num_eval_pts=0;
        //std::cout << "Number of eval points incoming: " << eval_pts_incoming.size() << std::endl;
        eval_points = boost::shared_array<EvalPt>(new EvalPt[eval_pts_incoming.size()]);
        assert(eval_points.get() != NULL); // check the allocation succeeded
        for (std::vector<Vector>::const_iterator it=eval_pts_incoming.begin(), end=eval_pts_incoming.end(); it != end; ++it)
        {
            eval_points[num_eval_pts].pt() = *it;
            eval_points[num_eval_pts].set_idx(num_eval_pts);
            ++num_eval_pts;
        }
        
        // Get the octree
        const ParallelTree& tree = *(ParallelFMMOctreeProxy.ckLocalBranch());
        typedef ParallelTree::NodeT NodeT;
        
        typedef std::map<OctreeIndexer, std::vector<EvalPt> > NodeEvalMap;
        NodeEvalMap mapper;
        
        // figure out which node each evaluation point lies in, then send 
        // work to the corresponding fmm_worker
        for (size_t ii=0; ii < num_eval_pts; ++ii)
        {
            const NodeT& node = tree.get_node(eval_points[ii]);
            NodeEvalMap::iterator find_it = mapper.find(node.get_idx());
            if (find_it == mapper.end())
            {
                mapper[node.get_idx()] = std::vector<EvalPt>();
                find_it = mapper.find(node.get_idx());
                assert(find_it != mapper.end());
            }
            find_it->second.push_back(EvalPt());
            *(find_it->second.rbegin()) = eval_points[ii];
        }
        
        // Create the messages to send work to the corresponding FMM node
        for (NodeEvalMap::const_iterator it=mapper.begin(), end=mapper.end(); it != end; ++it)
        {
            // Create a message
            const std::vector<EvalPt>& eval_list = it->second;
            
            int sizes[2];
            sizes[0] = eval_list.size();
            sizes[1] = 1;
            Eval_Message* msg = new(sizes, 0) Eval_Message;
            msg->length = eval_list.size();
            EvalPt* data = msg->data;
            for (size_t ii=0; ii < eval_list.size(); ++ii)
            {
                EvalPt* reinit_ep = new(&(data[ii])) EvalPt();
                data[ii].init();
                data[ii] = eval_list[ii];
            }
            
            msg->cb = CkCallback(CkIndex_Vanilla_FMM_Evals::receive_eval_results(NULL), thisProxy);
            msg->target = it->first;
            
            //std::cout << "Sending work to: " << it->first << std::endl;
            message_queue.push_back(msg);
            //VanillaFMMWorkerProxy[it->first].evaluate(msg);
            
            // increment the list of pending items
            ++pending;
        }
        
        // Now wait until somebody calls evaluate
        
    }
    
    void evaluate(CkCallback cb) 
    {
        // stash the callback- it's what we trigger when we've done all evaluations
        stashed_cb = cb;

        for (std::vector<Eval_Message*>::iterator it=message_queue.begin(), end=message_queue.end(); it != end; ++it)
        {
             VanillaFMMWorkerProxy[(*it)->target].evaluate(*it);
        }
        message_queue.clear();
        
        // now we chill out for a while and wait for the results to come back via
        // receieve_eval_results calls.  When the pending counter hits zero, all
        // the work is done and we send the results back in a message. (see below)
        

    }
    
    void receive_eval_results(Eval_Message* msg)
    {
        //std::cout << "Received some results... (pending=" << pending << ")" << std::endl;
        
        // each eval point contains a unique index corresponding to its position in the 
        // eval_points array
        for (size_t ii=0; ii < msg->length; ++ii)
        {
            // charm++ should guarantee only one entry method being called at a time
            // on any given chare, so shouldn't need thread safety here
            size_t idx = msg->data[ii].get_idx();
            eval_points[idx] += msg->data[ii]; // the += operator is not thread-safe
        }
        
        delete msg;
        
        // trigger the callback if we have received all the results
        if (--pending == 0)
        {
            // copy the eval points into a results message
            int sizes[2];
            sizes[0] = num_eval_pts;
            sizes[1] = 1;
            Eval_Message* results_msg = new(sizes, 0) Eval_Message;
            results_msg->length = num_eval_pts;
            EvalPt* data = results_msg->data;
            for (size_t ii=0; ii < num_eval_pts; ++ii)
            {
                EvalPt* reinit_ep = new(&(data[ii])) EvalPt();
                data[ii].init();
                data[ii] = eval_points[ii];
            }
            stashed_cb.send(results_msg);
            
            // minor optimisation: clear the eval_points memory
            eval_points.reset();
            
        }
    }
    
    Vanilla_FMM_Evals(CkMigrateMessage *msg) {
        std::cerr << "Haven't written a pup function for Vanilla_FMM_Evals!" << std::endl;
        throw std::exception();
    }
    
private:

    size_t pending;
    CkCallback stashed_cb;
    size_t num_eval_pts;
    boost::shared_array<EvalPt> eval_points;
    std::vector<Eval_Message*> message_queue;
    CProxy_VanillaFMMWorker VanillaFMMWorkerProxy;

};

} // end namespace

#endif // __VANILLA_FMM_EVALS_H
