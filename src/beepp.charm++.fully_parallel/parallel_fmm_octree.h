#ifndef __PARALLEL_FMM_OCTREE_H__
#define __PARALLEL_FMM_OCTREE_H__

#include "../common/charm_prereq.h"
#include <string>
#include <vector>
#include "../fmm/octree.h"

template<typename CType, typename CProxy_FMMWorkerT>
class ParallelFMMOctree : public CBase_ParallelFMMOctree<CType,CProxy_FMMWorkerT>, public Octree<Node<CType>,CType>
{

public:

    typedef CBase_ParallelFMMOctree<CType,CProxy_FMMWorkerT> super_charm;
    typedef Octree<Node<CType>,CType> super;

    /// Constructors ///
    ParallelFMMOctree(unsigned int max_items_per_node,
                        Vector universe_centre,
                        double universe_edge_length, 
                        unsigned int expected_total=0,
                        CkCallback cb_completed=CkCallback(CkCallback::ignore)) : super(max_items_per_node, universe_centre, universe_edge_length), expect(expected_total), stashed_cb(cb_completed)
    {
        global_lock = CmiCreateLock();
        if (expected_total == 0) {
            super_charm::contribute(cb_completed);
        }
    }
                        
    // should never migrate- it's a nodegroup object.
    ParallelFMMOctree(CkMigrateMessage *msg) : super(0, Vector(0,0,0), 0.0) { throw std::exception(); }

    ~ParallelFMMOctree() { 
        clear_waves();
        CmiDestroyLock(global_lock); 
    }

    /// Entry methods
    void insert(std::vector<CType> insertions) {
        CkGroupID gid = super_charm::thisgroup;
        //std::cout << "ParallelFMMOctree # " << gid.idx << " is inserting " << insertions.size() << " objects" << std::endl;
        for (typename std::vector<CType>::const_iterator it=insertions.begin(), end=insertions.end();
                it != end;
                ++it)
        {
            super::add(*it);
        }
        
        // check filledness
        if (super::size() == expect) {
            finalize(stashed_cb);
        }
    }
    
    void insert(CType single) {
         super::add(single);
        
         // check filledness
         if (super::size() == expect) {
             finalize(stashed_cb);
        }
    }
    
    void finalize(CkCallback cb)
    {
        super::build_neighbourhoods();
        //super::optimize_neighbourhoods();
        super::remove_empty_nodes();
        super_charm::contribute(cb);
    }
    
    template<typename T> // MUST be a Charm++ message type passed in- otherwise we'll leak memory (a *lot*)
    void get_waves_lock(const OctreeIndexer& idx, T* &out)
    {
        CmiLock(global_lock);
        bool done_lock=false;
        
        BufferMap::iterator find_it = buffers.find(idx);
        if (find_it == buffers.end())
        {
            T* new_ptr = new T(); // construct (i sure hope this is a message type)
            //std::cout << "Sizeof(T) is: " << sizeof(T) << std::endl;
            buffers[idx] = new_ptr;
            find_it = buffers.find(idx);
        }
        
        LockMap::iterator lock_it = locks.find(idx);
        if (lock_it == locks.end())
        {
            locks[idx] = CmiCreateLock();
            lock_it = locks.find(idx);
            CmiLock(lock_it->second);
            done_lock = true;
        }
        
        CmiUnlock(global_lock);
        
        // lock the buffer -- will be released by corresponding call to release_waves_lock
        if (!done_lock) {
            CmiLock(lock_it->second);
        }
        
        // Argh- nasty evil kludge!  Oh well.
        out = reinterpret_cast<T*>(find_it->second);
        
        return;
        
    }
    
    void release_waves_lock(const OctreeIndexer& idx)
    {
        CmiLock(global_lock);
        LockMap::iterator lock_it = locks.find(idx);
        assert(lock_it != locks.end()); // should always exist
        CmiUnlock(lock_it->second);
        CmiUnlock(global_lock);
        
        return;
    }
    
    void clear_waves() {
        
        CmiLock(global_lock);
        
        for (LockMap::iterator lock_it=locks.begin(), lock_end=locks.end(); lock_it != lock_end; ++lock_it)
        {
            // should not be locked
            assert( CmiTryLock(lock_it->second) == 0);
            CmiDestroyLock(lock_it->second);
        }
        locks.clear();
        
        for (BufferMap::iterator buf_it=buffers.begin(), buf_end=buffers.end(); buf_it != buf_end; ++buf_it)
        {
            delete buf_it->second;
        }
        buffers.clear();
        
        CmiUnlock(global_lock);
    }
    
    void request_data(CProxy_FMMWorkerT FMMWorkerProxy)
    {
        // need locks just in case anything is directly monkeying
        // about with the buffers - this makes this function
        // 'exclusive', guaranteed.
        CmiLock(global_lock);
        
        // get the entry
        for (BufferMap::iterator it=buffers.begin(), end=buffers.end(); it != end; ++it)
        {
            OctreeIndexer idx = it->first;
            
            LockMap::iterator lock_it = locks.find(idx);
            assert(lock_it != locks.end());
            assert(CmiTryLock(lock_it->second) == 0);
        
            // send the data
            if (FMMWorkerProxy[idx].ckLocal() != NULL)
            {
                FMMWorkerProxy[idx].ckLocal()->receive_incoming_wave(it->second); // yuk
            }
            else
            {
                FMMWorkerProxy[idx].receive_incoming_wave(it->second); // even more yuk
            }
            
            // destroy the lock
            CmiDestroyLock(lock_it->second);
        }
        
        buffers.clear();
        locks.clear();
        
        // release global lock
        CmiUnlock(global_lock);
        
        return;
        
    }
    
private:

    unsigned int expect;
    CkCallback stashed_cb;
    
    CmiNodeLock global_lock;
    
    typedef std::map< OctreeIndexer, CkMessage* > BufferMap;
    typedef std::map< OctreeIndexer, CmiNodeLock > LockMap;
    BufferMap buffers;
    LockMap locks;
    
};
#endif // __PARALLEL_FMM_OCTREE_H__
