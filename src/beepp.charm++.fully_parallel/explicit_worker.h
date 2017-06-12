#ifndef __EXPLICIT_WORKER_H__
#define __EXPLICIT_WORKER_H__

#include "prerequisites.h"
#include <boost/shared_array.hpp>
#include <boost/scoped_array.hpp>
#include "../bem/local_integrations.h"

class ExplicitWorker : public CBase_ExplicitWorker
{

public:

    /// Constructors ///
    ExplicitWorker();
    ExplicitWorker(CkMigrateMessage *msg);
    virtual ~ExplicitWorker() {}
    
    /// Entry methods
    void precalc(double kappa, unsigned int total_num_patches);
    void evaluate(double kappa);

    void pup(PUP::er &p);
    void triggerLB() {
      
#ifdef OPENCL 
#ifdef CACHE_GPU_RESULTS
      possibly_migrating = true;
      AtSync();
#endif
#endif
      
    }
    void ResumeFromSync() {
      possibly_migrating = false;
    }

private:

    void peak_integral(const CharmNodePatch& src_patch, const CharmNodePatch& targ, double& fpeak, double& hpeak, double kappa);

    std::vector<LintArray_Size> explicit_integrations;
    
#ifdef USE_JUFFER_PEAKS
    boost::scoped_array<double> juffer_peaks;
#endif
    bool possibly_migrating;
    bool done_precalcs;
};

#endif //__EXPLICIT_WORKER_H__
