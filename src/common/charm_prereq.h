#ifndef __CHARM_PREREQ_H
#define __CHARM_PREREQ_H

#include "../fmm/interaction_list.h"
#include "../common/matrix_types.h"
#include "../common/octree_indexer.h"
#include "../common/multipole_holder.h"
#include "../fmm/fmm_globals.h"
#include "../common/charge.h"
#include "../fmm/eval_pt.h"
#include <charm++.h>
#include <comlib.h>
#include <pup_stl.h>

using fmm::EvalPt;

static void resume_thread(void* ptr, double t)
{
    CthAwaken((CthThreadStruct*) ptr);
}

class CkArrayIndexOctreeIndexer : public CkArrayIndex {

public:

    CkArrayIndexOctreeIndexer()
    {
        OctreeIndexer *idx = new(index) OctreeIndexer();
        nInts=sizeof(*idx)/sizeof(int);
        assert(nInts == 3);
    }

    CkArrayIndexOctreeIndexer(const OctreeIndexer &in)
    {
        OctreeIndexer *idx = new(index) OctreeIndexer(in);
        nInts=sizeof(*idx)/sizeof(int);
        assert(nInts == 3);
    }

    // copy c'tor
    CkArrayIndexOctreeIndexer(const CkArrayIndexOctreeIndexer &other)
    {
        OctreeIndexer *idx = new(index) OctreeIndexer(other);
        nInts=sizeof(*idx)/sizeof(int);
        assert(nInts == 3);
    }

    operator OctreeIndexer &() {return *reinterpret_cast<OctreeIndexer*>(CkArrayIndex::index);}
    operator const OctreeIndexer &() const {return *reinterpret_cast<const OctreeIndexer*>(CkArrayIndex::index);}

    inline bool operator<(const OctreeIndexer& other) const
    {
        return static_cast<const OctreeIndexer&>(*this).as_hash_number() < other.as_hash_number();
    }

    /// PUP Routine ///
    inline void pup(PUP::er &p) {
        p(CkArrayIndex::index,3);
    }

};

class CharmEvalPt
{
    
public:
    
    CharmEvalPt() : pt(0,0,0), pot(0), field(0,0,0), idx(0) {}
    
    Vector pt;
    double pot;
    Vector field;
    size_t idx;
    
    /// PUP Routine ///
    inline void pup(PUP::er &p) {
        p | pt;
        p | pot;
        p | field;
        p | idx;
    }
};

#endif // __CHARM_PREREQ_H


