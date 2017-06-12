#ifndef __CHARM_GMRES_H__
#define __CHARM_GMRES_H__

#include "prerequisites.h"
#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include "rhs_handler.h"
#include <string>

class GMRES : public CBase_GMRES
{

public:

    /// Constructors ///
    GMRES();
    GMRES(CkMigrateMessage *msg);

    void solve(std::string output_filename,
				double kappa,
				double Dsolvent,
				unsigned int,
				FH_Values rhs);

    void pup(PUP::er &p) {
    	CBase_GMRES::pup(p);
	}
    
private:

};

#endif // __CHARM_GMRES_H__
