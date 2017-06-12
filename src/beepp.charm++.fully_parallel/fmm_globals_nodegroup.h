/*
 * fmm_globals_nodegroup.h
 *
 *  Created on: 6 Sep 2010
 *      Author: david
 */

#ifndef FMM_GLOBALS_NODEGROUP_H_
#define FMM_GLOBALS_NODEGROUP_H_

#include "prerequisites.h"
#include "charm++.h"
#include <vector>
#include <boost/shared_ptr.hpp>

#include "fmm_globals_nodegroup.decl.h"

template<int NTERMS, int NLAMBS>
class FMM_Globals_NodeGroupT : public CBase_FMM_Globals_NodeGroupT<NTERMS, NLAMBS>
{

public:

	 FMM_Globals_NodeGroupT() {
		 	 lock = CmiCreateLock();
	 }

	 ~FMM_Globals_NodeGroupT() {
	 	CmiDestroyLock(lock);
	 }

	 const fmm::Level_Dependent_FMM_Globals<NTERMS, NLAMBS>& get_globs_for_beta_level(double beta, int level, double universe_edge_length)
	 {
	 	CmiLock(lock);
	 	while (level >= level_globs_cache_beta.size())
	 	{
	 		int lev = level_globs_cache_beta.size();
	 		GlobPtr ptr(new GlobT(beta, lev, universe_edge_length));
	 		level_globs_cache_beta.push_back(ptr);
	 	}
	 	CmiUnlock(lock);
	 	return *(level_globs_cache_beta[level]);
	 }

	 const fmm::Level_Dependent_FMM_Globals<NTERMS, NLAMBS>& get_globs_for_beta0_level(double beta0, int level, double universe_edge_length)
	 {
	 	CmiLock(lock);
	 	while (level >= level_globs_cache_beta0.size())
	 	{
	 		int lev = level_globs_cache_beta0.size();
	 		GlobPtr ptr(new GlobT(beta0, lev, universe_edge_length));
	 		level_globs_cache_beta0.push_back(ptr);
	 	}
	 	CmiUnlock(lock);
	 	return *(level_globs_cache_beta0[level]);
	 }

	 inline const fmm::FMM_Globals<NTERMS>& get_fmm_globs() { return fmm_globs; };

private:

	 // NB: storing shared pointers to the FMM Globals struct is better than stashing the actual
	 // thing into the std::vector for 2 reasons:
	 // 1) it saves a complicated copy-construct (not written) of the Level_Dependent_FMM_Globals struct
	 // 2) it is thread safe in that references to the actual data obtained by the public member functions
	 // will always be valid, even if the std::vector concurrently gets resized to fit some more data in,
	 // because only the memory pointer of the shared_ptr will move, not the location pointed to by
	 // the shared pointer itself.
	 typedef fmm::Level_Dependent_FMM_Globals<NTERMS, NLAMBS> GlobT;
	 typedef boost::shared_ptr<GlobT> GlobPtr;
	 std::vector< GlobPtr > level_globs_cache_beta;
	 std::vector< GlobPtr > level_globs_cache_beta0;
	 fmm::FMM_Globals<NTERMS> fmm_globs;

	 CmiNodeLock lock;


};

#define CK_TEMPLATES_ONLY
#include "fmm_globals_nodegroup.def.h"
#undef CK_TEMPLATES_ONLY

#endif /* FMM_GLOBALS_NODEGROUP_H_ */
