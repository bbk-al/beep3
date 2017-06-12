#ifndef __VANILLA_FMM_H
#define __VANILLA_FMM_H

#include "opencl_nodegroup.decl.h"
#include "opencl_nodegroup.h"
#include "../common/charm_prereq.h"
#include <string>
#include <vector>
#include <algorithm>
#include "../fmm/fmm_octree.h"
#include "../fmm/eval_pt.h"

extern /* readonly */ CProxy_OpenCL_NodeGroup OpenCL_NodeGroupProxy;

class EvalMessage : public CMessage_EvalMessage {

public:

    int length;
    int length_cb;
    int length_user_flag;
    EvalPt* data;
    CkCallback cb;
    int user_flag;

};

template<int NTERMS, int NLAMBS, int NWAVES> 
         class VanillaFMM : public CBase_VanillaFMM<NTERMS, NLAMBS, NWAVES>, 
                                   public fmm::FMM_Octree_T<Charge, NTERMS, NLAMBS, NWAVES>
{

public:

    typedef CBase_VanillaFMM<NTERMS, NLAMBS, NWAVES> super_charm;
    typedef fmm::FMM_Octree_T<Charge,NTERMS,NLAMBS,NWAVES> super;
    typedef typename super::FMM_Node NodeT;

    /// Constructors ///
    VanillaFMM(unsigned int max_items_per_node,
                        Vector universe_centre,
                        double universe_edge_length, 
                        unsigned int expected_total=0,
                        CkCallback cb_completed=CkCallback(CkCallback::ignore)) : super(max_items_per_node, universe_centre, universe_edge_length), expect(expected_total), stashed_cb(cb_completed)
    {
        if (expected_total == 0) {
            super_charm::contribute(cb_completed);
        }
        my_ctr = 0;
        kappa = -1;
    }

    // should never migrate- it's a nodegroup object.
    VanillaFMM(CkMigrateMessage *msg) : super(0, Vector(0,0,0), 0.0) { throw std::exception(); }

    ~VanillaFMM() { }

    /// Entry methods
    void insert(std::vector<Charge> insertions) {
        CkGroupID gid = super_charm::thisgroup;
        //std::cout << "VanillaFMM # " << gid.idx << " is inserting " << insertions.size() << " objects" << std::endl;
        for (typename std::vector<Charge>::const_iterator it=insertions.begin(), end=insertions.end();
                it != end;
                ++it)
        {
            super::add(*it);
        }
        
        // check filledness
        if (super::size() == expect) {
            super_charm::thisProxy[CkMyPe()].finalize(stashed_cb);
        }
    }
    
    void insert(Charge single) {
         super::add(single);
        
         // check filledness
         if (super::size() == expect) {
             super_charm::thisProxy[CkMyPe()].finalize(stashed_cb);
         }
    }
  
    void callback_when_filled(unsigned int num, CkCallback cb) {
         
        // send callback when reach filledness
        expect = num;
        stashed_cb = cb;
        if (super::size() >= expect) {
             super_charm::thisProxy[CkMyPe()].finalize(stashed_cb);
        }
        
    }
    
    void finalize(CkCallback cb)
    {
        // delete nodes which this compute group is not interested in
        decision_level = level_with_more_nodes_than_procs();
        //super::linearize(decision_level);
        int total_desc=0;
        typename super::NodeList& nlist = super::get_node_list(decision_level);
        
        int erased=0;
#if 1
        int num_pts_per_pe = static_cast<int>(floor(static_cast<double>(super::get_everything().size()) / CkNumPes()));
        int pe_ctr=0, nctr=0;
        for (typename super::NodeList::iterator it=nlist.begin(), end=nlist.end(); it != end; ++it)
        {
            // add number of points
            int new_size = nctr + it->second->size(); 
            if (new_size > num_pts_per_pe && (new_size - num_pts_per_pe) > (num_pts_per_pe - nctr))
            {
                ++pe_ctr;
                if (pe_ctr >= CkNumPes()) { 
                    pe_ctr = CkNumPes()-1; 
                    nctr += it->second->size();
                }
                else
                {
                    nctr = it->second->size();
                }
                
            }
            else
            {
            	nctr += it->second->size();
            }
            
            // store which processor this node is evaluated on
            NodeToProc[it->second->get_idx()] = pe_ctr;

            if (pe_ctr != CkMyPe())
            {
                erased += super::erase_children(*(it->second));
                it->second->unset_is_leaf();
                it->second->set_shadow();
            }
        }
        
         //std::cout << "Total erased: " << erased << std::endl;
#else
        
        // get nodes in size order
        std::vector<NodeT*> nodes;
        for (typename super::NodeList::iterator it=nlist.begin(), end=nlist.end(); it != end; ++it)
        {
            nodes.push_back(it->second);
        }
        std::sort(nodes.begin(), nodes.end(), typename super::NodeSorter());
        
        int pe_ctr=0;
        for (typename std::vector<NodeT*>::iterator nit=nodes.begin(), nend=nodes.end(); nit != nend; ++nit)
        {
            // store which processor this node is evaluated on
            NodeToProc[(*nit)->get_idx()] = pe_ctr;
            
            if (pe_ctr != CkMyPe())
            {
                super::erase_children(**nit);
                (*nit)->unset_is_leaf();
                (*nit)->set_shadow();
            }
            ++pe_ctr;
            if (pe_ctr >= CkNumPes()) { pe_ctr = 0; }
        }
#endif

        super::build_neighbourhoods();
        super::remove_empty_nodes();
        
        super_charm::contribute(cb);
    }
    
    void solve(double kappa_in, CkCallback cb)
    {
        kappa = kappa_in;
        super::solve(kappa);
        super_charm::contribute(cb);
    }

    void solve_and_evaluate(double kappa_in, CkCallback cb)
    {
        kappa = kappa_in;
        super::solve(kappa);
        //std::cout << CkMyPe() << " FMM timings: " << super::get_timing_info() << std::endl;
        super_charm::thisProxy[CkMyPe()].evaluate(cb);
    }

    unsigned short level_with_more_nodes_than_procs() const
    {
        for (unsigned short level=super::get_top_level(); level <= super::get_bottom_level(); ++level)
        {
            if (super::get_num_nodes_on_level(level) >= CkNumPes()) 
            {
                return level;
            }
            
        }
        return super::get_bottom_level();
    }
    
    void evaluate(CkCallback cb_results)
    {
        if (kappa == -1) { 
            std::cerr << "Warning: kappa not set in VanillaFMM::evaluate.  Have you called solve?" << std::endl;
        }
        
        std::vector<EvalPt*> pts;
        typedef typename NodeT::ContentList All;
        All& everything = super::get_everything();
        size_t ctr=0;
        for (typename All::const_iterator it=everything.begin(), end=everything.end(); it != end; ++it)
        {
            const NodeT& node = super::get_node(**it, decision_level);
            OctreeIndexer idxer = node.get_idx();
            unsigned short eval_level = idxer.get_level();
            assert(eval_level <= decision_level);
            if (  (eval_level < decision_level && (idxer.as_hash_number() % CkNumPes() == CkMyPe())) 
                  || NodeToProc[idxer] == CkMyPe() ) 
            {
                // any processor can handle it if the point is higher than the level
                // at which we have split work across nodes. Choose processor
                // based on the hash number (only one processor actually does the work)
                pts.push_back(new EvalPt(static_cast<Vector&>(**it)));
                (**pts.rbegin()).set_idx(ctr);
            }
            ++ctr;
        }
        
#ifdef OPENCL
        OpenCL_Handler& ocl_handler = OpenCL_NodeGroupProxy.ckLocalBranch()->get_ocl_handler();
        super::evaluate_many(pts, ocl_handler);
#else
        super::evaluate_many(pts);
#endif
        
        // delete allocated memory
        for (std::vector<EvalPt*>::iterator it=pts.begin(), end=pts.end(); it != end; ++it)
        {
            if ((**it).get_idx() <= 10) {
                std::cout << **it << std::endl;
            }
            delete *it;
        }
        pts.clear();
        
        super_charm::contribute(cb_results);
        
        return;
    }
    
    // this should be sent only to the local chare
    void __private_evaluate(EvalMessage *msg)
    {
        assert(kappa != -1);
        int num_evals = msg->length;
        EvalPt* incoming = msg->data;
        std::vector<EvalPt*> pts;
        for (int ii=0; ii < num_evals; ++ii)
        {
            pts.push_back(&(incoming[ii]));
        }

#ifdef OPENCL
        OpenCL_Handler& ocl_handler = OpenCL_NodeGroupProxy.ckLocalBranch()->get_ocl_handler();
        super::evaluate_many(pts, ocl_handler);
#else
        super::evaluate_many(pts);
#endif

        // trigger callback
        msg->cb.send(msg);
    }
    
    // this should be sent only to the local chare
    void evaluate(EvalMessage *msg)
    {
        int num_evals = msg->length;
        EvalPt* incoming = msg->data;
        
        std::vector<EvalPt*> pts;
        std::map<int, std::vector<int> > outbound;
        
        for (int ctr=0; ctr < num_evals; ++ctr)
        {
            EvalPt* ep = &(incoming[ctr]);
            
            const NodeT& node = super::get_node(*ep, decision_level);
            OctreeIndexer idxer = node.get_idx();
            unsigned short eval_level = idxer.get_level();
            assert(eval_level <= decision_level);
            //std::cout << ep->pt() << " @ " << idxer << " " << eval_level << " " << decision_level << std::endl;
            int proc = (eval_level < decision_level) ? idxer.as_hash_number() % CkNumPes() : NodeToProc[idxer];
            if (outbound.find(proc) == outbound.end()) {
                outbound[proc] = std::vector<int>();
            }
            outbound[proc].push_back(ctr);
        }
        
        // send outbound messages
        for (std::map<int, std::vector<int> >::const_iterator it=outbound.begin(), end=outbound.end();
            it != end; ++it)
        {
            int target_proc = it->first;
            const std::vector<int>& idxs = it->second;
            EvalMessage* msg_out = new (idxs.size(), 0) EvalMessage;
            msg_out->length = idxs.size();
            
            int ctr=0;
            for (std::vector<int>::const_iterator idx_it=idxs.begin(), idx_end=idxs.end(); idx_it != idx_end; ++idx_it)
            {
                EvalPt& ep = incoming[*idx_it];
                EvalPt* dummy = new (&(msg_out->data[ctr++])) EvalPt(ep);
                dummy->set_idx(*idx_it);
            }
            
            msg_out->user_flag = my_ctr;
            msg_out->cb = CkCallback(CkIndex_VanillaFMM<18,18,300>::recv_data(NULL), CkMyPe(), super_charm::thisProxy);
            super_charm::thisProxy[target_proc].__private_evaluate(msg_out);
        }
        
        // this keeps track of the work as it gets done    
        pending[my_ctr++] = Pending(msg, num_evals);
        
        return;
    }
    
    void recv_data(EvalMessage *msg)
    {
        int n = msg->user_flag;
        assert(pending.find(n) != pending.end());
        Pending& p = pending[n];
        EvalMessage* return_msg = p.first;
        int& evals_pending = p.second;
        
        // process the incoming message - store evals in the
        // stored message (stored in pending list)
        for (int ii=0; ii < msg->length; ++ii)
        {
            EvalPt& ep = msg->data[ii];
            return_msg->data[ep.get_idx()] += ep;
        }
        
        // if we have now completed the pending evaluations
        // we can return the message to the original sender
        evals_pending -= msg->length;
        if (evals_pending == 0)
        {
            return_msg->cb.send(return_msg);
            pending.erase(n);
        }
        
        delete msg;
        
    }
    
private:

    unsigned int expect;
    CkCallback stashed_cb;
    
    typedef std::pair<EvalMessage*, int> Pending;
    std::map<int, Pending> pending;
    
    double kappa;
    unsigned short decision_level;
    std::map<OctreeIndexer, int> NodeToProc;
    int my_ctr;
    
};

#endif // __VANILLA_FMM_H
