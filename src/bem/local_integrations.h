/*
* local_integrations.h
*
*  Created on: 30 Jul 2010
*      Author: david
*/

#ifndef LOCAL_INTEGRATIONS_H_
#define LOCAL_INTEGRATIONS_H_

#include "bem_kernels.h"

#ifdef OPENCL
#include <CL/opencl.h>
#include "opencl_bem_structs.h"
#else
typedef struct ocl_lint
{
    unsigned int src_global_idx;
    unsigned int targ_global_idx;
    float Apt;
    float Bpt;
    float Cpt;
    float Dpt;
} OpenCL_LocalIntegrationResult;
#endif

#include "node_patch.h"
#include "boost/shared_array.hpp"
#ifdef __CHARMC__
#include <charm++.h>
#endif

class LocalIntegrations : public OpenCL_LocalIntegrationResult
{

public:

    typedef OpenCL_LocalIntegrationResult super;

    LocalIntegrations() {
        super::src_global_idx = 0;
        super::targ_global_idx = 0;
        super::Apt = 0;
        super::Bpt = 0;
        super::Cpt = 0;
        super::Dpt = 0;
    }
    
    LocalIntegrations(double kappa,
                        const BasicNodePatch& src_patch,
                        const BasicNodePatch& targ_patch)
    {
        init(kappa, src_patch, targ_patch);
    }
    
    LocalIntegrations(unsigned int src_idx,
                        unsigned int targ_idx,
                        float A,
                        float B,
                        float C,
                        float D)
    {
        super::src_global_idx = src_idx;
        super::targ_global_idx = targ_idx;
        super::Apt = A;
        super::Bpt = B;
        super::Cpt = C;
        super::Dpt = D;
    }    

    // copy constructor
    LocalIntegrations(const LocalIntegrations& other) {

        super::src_global_idx = other.src_global_idx;
        super::targ_global_idx = other.targ_global_idx;
        super::Apt = other.Apt;
        super::Bpt = other.Bpt;
        super::Cpt = other.Cpt;
        super::Dpt = other.Dpt;
    }

    inline void set(unsigned int src_idx,
            unsigned int targ_idx,
            float A,
            float B,
            float C,
            float D)
    {
        super::src_global_idx = src_idx;
        super::targ_global_idx = targ_idx;
        super::Apt = A;
        super::Bpt = B;
        super::Cpt = C;
        super::Dpt = D;
    }   

    inline void init(double kappa,
              const BasicNodePatch& src_patch,
              const BasicNodePatch& targ_patch)
    {
        //const float epsilon = 1.0 / inv_epsilon;

        super::src_global_idx = src_patch.get_idx();
        super::targ_global_idx = targ_patch.get_idx();

        evaluate_local_BEM_kernels(static_cast<float>(kappa), src_patch, targ_patch, super::Apt, super::Bpt, super::Cpt, super::Dpt);

    }

    inline bool insane(double max_value)
    {
        if (fabs(super::Apt) > max_value || 
            fabs(super::Bpt) > max_value || 
            fabs(super::Cpt) > max_value || 
            fabs(super::Dpt) > max_value)
        {
            return true;
        }
        return false;
    }

    inline void evaluate_local_contributions(const double fvals[], const double hvals[], double fresults[], double hresults[]) const
    {
        // multiply the pre-calc'd BEM kernel matrix elements by the correct bits of the lhs vector
        fresults[super::src_global_idx] += super::Bpt * fvals[super::targ_global_idx] - super::Apt * hvals[super::targ_global_idx];
        hresults[super::src_global_idx] += super::Dpt * fvals[super::targ_global_idx] - super::Cpt * hvals[super::targ_global_idx];

        return;
    }

    // version using kahan compensated addition
    inline void evaluate_local_contributions(const double fvals[], const double hvals[],
                                                double fresults[], double hresults[],
                                                double kahan_f[], double kahan_h[]) const
    {
        // multiply the pre-calc'd BEM kernel matrix elements by the correct bits of the lhs vector
        {
            double& sum = fresults[super::src_global_idx];
            double val = super::Bpt * fvals[super::targ_global_idx] - super::Apt * hvals[super::targ_global_idx];
            double& kahan = kahan_f[super::src_global_idx];
            {
                // kahan summation (aka compensated addition)
                double tmp_y = val - kahan;
                double tmp_t = sum + tmp_y;
                kahan = (tmp_t - sum) - tmp_y;
                sum = tmp_t;
            }
        }

        {
            double& sum = hresults[super::src_global_idx];
            double val = super::Dpt * fvals[super::targ_global_idx] - super::Cpt * hvals[super::targ_global_idx];
            double& kahan = kahan_h[super::src_global_idx];
            {
                // kahan summation (aka compensated addition)
                double tmp_y = val - kahan;
                double tmp_t = sum + tmp_y;
                kahan = (tmp_t - sum) - tmp_y;
                sum = tmp_t;
            }
        }

        return;
    }

    std::string str() const
    {
        std::ostringstream buf;
        buf << "From " << super::src_global_idx << " to " << super::targ_global_idx
                << " Apt=" << super::Apt << " Bpt=" << super::Bpt << " Cpt=" << super::Cpt << " Dpt=" << super::Dpt << " ";
        return buf.str();
    }

#ifdef __CHARMC__

    void pup(PUP::er &p)
    {
        p | super::src_global_idx;
        p | super::targ_global_idx;
        p | super::Apt;
        p | super::Bpt;
        p | super::Cpt;
        p | super::Dpt;
    }

#endif

};

inline std::ostream& operator<<(std::ostream& os, const LocalIntegrations& loc)
{
    os << loc.str();
    return os;
}

typedef boost::shared_array<LocalIntegrations> LintArray;
typedef std::pair< LintArray, size_t> LintArray_Size;

#endif /* LOCAL_INTEGRATIONS_H_ */
