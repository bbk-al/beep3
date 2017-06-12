/*
 * fh_values_nodegroup.cpp
 *
 *  Created on: 6 Sep 2010
 *      Author: david
 */

#include "prerequisites.h"
#include "fh_values_nodegroup.decl.h"
#include "fh_values_nodegroup.h"

#include <boost/shared_array.hpp>

FH_Values_NodeGroup::FH_Values_NodeGroup(unsigned int num_patches_) : num_patches(num_patches_), size(num_patches_*2)
{
	//std::cout << "Allocating room for " << size << " doubles" << std::endl;
	lhs = boost::shared_array<double>(new double[size]);
	results = boost::shared_array<double>(new double[size]);
	kahan = boost::shared_array<double>(new double[size]);
	assert(lhs.get() != NULL);
	assert(results.get() != NULL);
	assert(kahan.get() != NULL);

	lock = CmiCreateLock();
}

void FH_Values_NodeGroup::set(FH_Values vals, CkCallback cb)
{
	assert(lhs.get() != NULL);
	assert(vals.size() == size);
	assert(vals.get() != NULL);
	memcpy(lhs.get(), vals.get(), size*sizeof(double));
	memset(results.get(), 0, size*sizeof(double));
	memset(kahan.get(), 0, size*sizeof(double));
	contribute(cb);
}

#include "fh_values_nodegroup.def.h"
