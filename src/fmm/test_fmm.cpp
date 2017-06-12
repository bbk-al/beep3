// Test Program for FMM (C++ port by David Fallaize, based on Jingfang Huang's Fortan version, see fastmultipole.org)

#ifdef OPENCL
//#undef OPENCL
class OpenCL_Handler;
#endif

#include "fmm_octree.h"
#include "../common/math_vector.h"
#include "../common/matrix_types.h"
#include "../common/charge.h"
#include "eval_pt.h"
#include <vector>
#include <iostream>
#include <limits>
#include <cstdlib> 
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <sys/times.h>
#include "../common/useful_clock.h"

typedef std::numeric_limits< double > dbl;

const unsigned int default_num_repeats = 1;
const unsigned int default_max_objects = 100;
const double kappa = 0.1;
const unsigned int NUM_PTS_TO_COMPARE = 100;

static const double ticks_per_millisecond = sysconf(_SC_CLK_TCK) / 1000.;
static double subtract_time_elapsed(const tms& start, const tms& end)
{
    return ((end.tms_utime - start.tms_utime) + (end.tms_stime - start.tms_stime)) / ticks_per_millisecond;
}

int get_random_number(int low, int high)
{
    int range=(high-low)+1;
    return low+int(range*double(rand())/(RAND_MAX + 1.0));
}

int random_sign()
{
    return (rand() >= RAND_MAX/2) ? -1 : 1;
}

int main(int argc, char** argv)
{

    if (argc < 2) {
        std::cout << "Must specify a file containing charge locations in x,y,z,q\\n format" << std::endl;
        return 1;
    }
#ifdef OPENCL
    bool use_opencl = false;
    OpenCL_Handler global_ocl_handler;
#endif

    std::string filename(argv[1]);
    std::cout << "Hello.  This is the FMM test program.  Reading test data from " << filename << std::endl;
    
    unsigned int max_objects_per_cell = default_max_objects;
    if (argc >= 3) {
        max_objects_per_cell = boost::lexical_cast<unsigned int>(argv[2]);
        std::cout << "Using user-specified max_objects_per_cell=" << max_objects_per_cell << std::endl;
    
    }

    unsigned int num_repeats = default_num_repeats;
    if (argc >= 4) {
        num_repeats = boost::lexical_cast<unsigned int>(argv[3]);
    }

    srand((unsigned)time(0));
    std::cout.precision(dbl::digits10);
    
    std::ifstream fin;
    fin.open(filename.c_str());
    
    // check that the file stream is valid
    assert(fin.good());

    // list of charges
    std::vector<Charge> charges;

    // create list of evaluation points too
    std::vector<fmm::EvalPt_2ndDerivs*> eval_pts;

    Vector max(-1e99,-1e99,-1e99);
    Vector min(+1e99,+1e99,+1e99);
    double charge;
    Vector xyz;
    while(fin >> xyz.x >> xyz.y >> xyz.z >> charge)
    {
        // Now remove the extra stuff on the line you do not want.
        fin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
        
        max.x = (xyz.x > max.x) ? xyz.x : max.x;
        max.y = (xyz.y > max.y) ? xyz.y : max.y;
        max.z = (xyz.z > max.z) ? xyz.z : max.z;
        min.x = (xyz.x < min.x) ? xyz.x : min.x;
        min.y = (xyz.y < min.y) ? xyz.y : min.y;
        min.z = (xyz.z < min.z) ? xyz.z : min.z;

        charges.push_back(Charge(xyz,charge));
        eval_pts.push_back(new fmm::EvalPt_2ndDerivs(xyz) );
    }

    // close file
    fin.close();

    const Vector centre = (max + min) / 2.0;
    const double edge_length = (max - min).length() * 1.1;
    const unsigned int num_eval_pts = charges.size();
    std::cout << "Num Charges: " << charges.size() << " Centre: " << centre << " Edge length: " << edge_length << std::endl;

    double inserting = 0;
    double solving = 0;
    double fmm_eval = 0;
    double time_explicit_each = 0;
    double time_fmm = 0;
    double time_explicit = 0;
    double total_real = 0;
    double mean_err = 0;

    for (unsigned int repeat_ctr=0; repeat_ctr < num_repeats; ++repeat_ctr)
    {

    tms times_start_adding;
    tms times_done_adding;
    tms times_start_solving;
    tms times_done_solving;
    tms times_start_eval_cpu;
    tms times_done_eval_cpu;
    long time_start_eval_real;
    double time_eval_real;
    tms times_start_explicit;
    tms times_done_explicit;
    long total_start_real=myclock();

    fmm::FMM_Octree_6FIG_ACCURACY fmm(max_objects_per_cell, centre, edge_length);

    // hack max depth to measure quadratic scalings (non-FMM)
    //std::cout << "Hacking max depth to zero..." << std::endl;
    //fmm.set_max_depth(0);

    std::cout << "Adding charges..." << std::flush;
    times(&times_start_adding);
    int ctr=0;
    fmm.add(charges);
//    for (std::vector<Charge>::const_iterator it=charges.begin(), end=charges.end(); it != end; ++it)
//    {
//        fmm.add(*it);
//    }
    times(&times_done_adding);
    std::cout << "done" << std::endl;

    std::cout << "Solving FMM, kappa=" << kappa << std::endl;
    times(&times_start_solving);
    fmm.solve(kappa);
    times(&times_done_solving);

    std::cout << "Evaluating FMM at charge locations" << std::endl;
    time_start_eval_real = myclock();
    times(&times_start_eval_cpu);
#ifdef OPENCL
    if (use_opencl)
    {
        fmm.evaluate_many(eval_pts, global_ocl_handler);
    }
    else
    {
        fmm.evaluate_many(eval_pts);
    }
#else
    fmm.evaluate_many(eval_pts);
#endif
    times(&times_done_eval_cpu);
    time_eval_real = (myclock() - time_start_eval_real)/1e3;
    total_real += (myclock() - total_start_real)/1e3;

    unsigned int num_explicit = (NUM_PTS_TO_COMPARE > eval_pts.size()) ? eval_pts.size() : NUM_PTS_TO_COMPARE;

    if (repeat_ctr == 0)
    {
    // explicitly calculate the first 100 evaluation points
    std::cout << "Now doing it the slow O(N^2) way..." << std::endl;
    times(&times_start_explicit);
    double numer=0;
    double denom=0;
    for (std::vector<fmm::EvalPt_2ndDerivs*>::iterator ep_it=eval_pts.begin(), ep_end=ep_it+num_explicit; ep_it!=ep_end; ++ep_it)
    {
        const fmm::EvalPt_2ndDerivs& ep = **ep_it;
        fmm::EvalPt_2ndDerivs explicit_ep(static_cast<Vector&>(**ep_it));

        for (std::vector<Charge>::const_iterator ch_it=charges.begin(), ch_end=charges.end(); ch_it != ch_end; ++ch_it)
        {
            fmm::EvalPt_2ndDerivs::add_explicit_contrib(explicit_ep, *ch_it, kappa);
        
        }
        double explicit_pot = explicit_ep.get_potential();
        double e = fabs(explicit_pot - ep.get_potential());
        numer += e*e;
        denom += fabs(explicit_pot)*fabs(explicit_pot);
        
        std::cout << "FMM: " << ep << "\n" 
        		  << "Exp: " << explicit_ep << "\n\n";

    }
    times(&times_done_explicit);
    double err = sqrt(numer / denom);
    std::cout << "Error in potential is: " << err << std::endl;
    mean_err = err;
    time_explicit_each = subtract_time_elapsed(times_start_explicit, times_done_explicit) / static_cast<double>(num_explicit);

    }

    for (std::vector<fmm::EvalPt_2ndDerivs*>::iterator ep_it=eval_pts.begin(), ep_end=eval_pts.end(); ep_it!=ep_end; ++ep_it)
    {
        (**ep_it).reset();
    }

    inserting += subtract_time_elapsed(times_start_adding, times_done_adding);
    solving += subtract_time_elapsed(times_start_solving, times_done_solving);
    double t_explicit = fmm.get_timing_info().explicit_evaluations / 1000.;
    time_fmm += (time_eval_real - t_explicit);
    time_explicit += t_explicit;

    }

    std::cout << "Done." << std::endl;

    // timings
    inserting /= num_repeats;
    solving /= num_repeats;
    //fmm_eval /= num_repeats;
    time_fmm /= num_repeats;
    time_explicit /= num_repeats;
    total_real /= num_repeats;

    std::cout << "RESULTS: " << max_objects_per_cell << ", " << num_eval_pts << ", " << inserting << ", " << solving << ", " << time_fmm << ", " << time_explicit << ", " << time_explicit_each*num_eval_pts << ", " << total_real << ", " << mean_err << std::endl;

    // delete allocated memory for evaluation points
    std::cout << "Deleting allocated memory" << std::endl;
    for (std::vector<fmm::EvalPt_2ndDerivs*>::iterator it=eval_pts.begin(), end=eval_pts.end(); it != end; ++it)
    {
        delete *it;
    }
    eval_pts.clear();

    return 0;

}
