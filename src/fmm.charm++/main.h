#ifndef __MAIN_H__
#define __MAIN_H__
#include <boost/shared_ptr.hpp>

class Main : public CBase_Main
{

private:

    long start_time;
    double universe_edge_length;
    double kappa;

public:

    /// Constructors ///
    Main(CkArgMsg* msg);
    Main(CkMigrateMessage* msg);

    /// Entry Methods ///
    void create_workers();
    void fmm_worker_complete();

    // This gets called when all goes quiet
    void quiescenceHandler();

    // gets called on every processor at start
    static void initManualLB() {
        std::cout << "Disabling load balancing on processor " << CkMyPe() << std::endl;
        TurnManualLBOn();
    }

};

#endif //__MAIN_H__
