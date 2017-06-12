#ifndef __PREREQUISITES_H
#define __PREREQUISITES_H

// some commonly needed charm++ prerequites
// (commonly needed by me that is)
#include "../common/charm_prereq.h"

// we need the vanilla fmm prereqs in several places
#include "vanilla_fmm_prerequisites.h"

// the rest of this is specific to beepp
#include "../bem/mesh.h"
#include "../bem/mesh_instance.h"
#include "../bem/node_patch.h"
#include "../bem/local_integrations.h"

// Useful BEEP stuff like RunInfo and the default BEM neighbourhood size
#include "../beep/beep.h"

#include "charm_node_patch.h"

// A holder class for intermediate bem/fmm f,h values
class FH_Values
{

public:

	FH_Values() : _size(0), data(NULL) {}
	FH_Values(FH_Values& other) : _size(other._size), data(other.data) {}
	FH_Values(size_t size_, double* incoming) : _size(size_) {
		data = boost::shared_array<double>(new double[_size]);
		memcpy(data.get(), incoming, sizeof(double)*_size);
	}
	inline size_t size() const { return _size; }
	inline double* get() { return data.get(); }

	void pup(PUP::er &p) {

		p | _size;

		if (p.isUnpacking()) {
			data = boost::shared_array<double>(new double[_size]);
		}

        for (size_t ii=0; ii < _size; ++ii)
        {
            p | data[ii];
        }
    }

protected:
	size_t _size;
	boost::shared_array<double> data;

};

class EvalMessage;
template<int,int,int> class CProxy_VanillaFMM;
template<int,int,int,typename> class CProxy_BEM_FMM;

// This is jolly ugly i admit - these are the NTERM, NLAMBS, NWAVES
// template parms for the BEM-FMM class - 9/9/67 is 3 digit accuracy
#define NT 9
#define NL 9
#define NW 67

#endif
