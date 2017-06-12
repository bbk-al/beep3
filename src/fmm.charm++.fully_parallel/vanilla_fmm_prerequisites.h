#ifndef __VANILLA_FMM_PREREQUISITES_H
#define __VANILLA_FMM_PREREQUISITES_H

// some useful charm++ stuff and common includes
// Also the CkArrayIndexer analogue of the OctreeIndexer
// is in charm_prereq.h
#include "../common/charm_prereq.h"

// the rest of this is specific to the pfmm parallel fmm (vanilla)

// forward decls.
template<int , int > class CProxy_FMM_Globals_NodeGroupT;
template<typename CType, typename CProxy_FMMWorkerT> class CProxy_ParallelFMMOctree;

namespace vanilla_fmm
{

    template<int NTERMS> class MultipoleHolderT : public MultiHolder<1, BaseMultipoleHolder<NTERMS> > {};
    template<int NLAMBS, int NWAVES> class PlaneWaveHolderT : public MultiHolder<1, BasePlaneWaveHolder<NLAMBS, NWAVES> > {};
    
    // fwd declarations for some Proxy types so that we can define some useful typedefs
    template<int , int , int > class CProxy_VanillaFMMWorkerT;

    typedef CProxy_VanillaFMMWorkerT<18,18,300> CProxy_VanillaFMMWorker;
    typedef ::CProxy_FMM_Globals_NodeGroupT<18,18> CProxy_FMM_Globals_NodeGroup;
    typedef CProxy_ParallelFMMOctree<Charge, CProxy_VanillaFMMWorker> CProxy_ParallelTree;

}

#endif // __VANILLA_FMM_PREREQUISITES_H
