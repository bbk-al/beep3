/*
 * fh_values_nodegroup.h
 *
 *  Created on: 6 Sep 2010
 *      Author: david
 */

#ifndef FH_VALUES_NODEGROUP_H_
#define FH_VALUES_NODEGROUP_H_

#include "prerequisites.h"
#include <vector>
#include <boost/shared_array.hpp>

class FH_Values_NodeGroup : public CBase_FH_Values_NodeGroup {

public:

	FH_Values_NodeGroup(unsigned int);
	~FH_Values_NodeGroup() {
		CmiDestroyLock(lock);
	}

	inline const double* fvals() const { assert(lhs.get() != NULL); return &(lhs[0]); }
	inline const double* hvals() const { assert(lhs.get() != NULL); return &(lhs[num_patches]); }

	inline const double* get_results() const { assert(results.get() != NULL); return results.get(); }

	void add_results(double* add_me)
	{

		// Get lock
		CmiLock(lock);

		for (unsigned int ctr=0; ctr < size; ++ctr)
		{
			// kahan summation (aka compensated addition- kahan array holds low-order part)
			double tmp_y = add_me[ctr] - kahan[ctr];
			double tmp_t = results[ctr] + tmp_y;
			kahan[ctr] = (tmp_t - results[ctr]) - tmp_y;
			results[ctr] = tmp_t;
		}

		// Unlock
		CmiUnlock(lock);
	}

	void set(FH_Values vals, CkCallback cb);

	void reduce(CkCallback cb)
	{
		// Reduce elementwise down the entire array -- use cb as callback (should point
		// to done_results_reduction in Main chare).
		contribute(size*sizeof(double), results.get(), CkReduction::sum_double, cb);
	}

	inline unsigned int get_num_patches() const { return num_patches; }

private:

	// for thread safety when adding results
	CmiNodeLock lock;

	const unsigned int num_patches;
	const size_t size;
	boost::shared_array<double> lhs;
	boost::shared_array<double> results;
	boost::shared_array<double> kahan;

};

#endif /* FH_VALUES_NODEGROUP_H_ */
