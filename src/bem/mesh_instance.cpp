/*      Author: david fallaize    Created on: 27 Jul 2010 */
/*      Modified: adam light    on: 9 Mar 2012  */

/*! \file mesh_instance.cpp
 * \brief This module implements the mesh instance and list classes.
 *
 */

#include "mesh_instance.h"
#include "constants.h"
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <functional>
#include <iostream>
#include <iomanip>
#include <unordered_set>
#include <random>
#include <boost/functional/hash.hpp>


MeshInstance::MeshInstance(
	unsigned int mesh_lib_id,
	unsigned int mesh_instance_id,
	MeshList& mesh_library,
	const Vector& offset,
	const Quaternion& rot,
	double protein_dielectric,
	double solvent_dielectric,
	unsigned int num_quad_points,
	unsigned int num_qual_points,
	bool _silent
	) : instance_id(mesh_instance_id),
		silent(_silent),
		xyz_offset(offset),
		rotation(rot),
		psel(0),
		quad_points_per_triangle(num_quad_points),
		qual_points_per_triangle(num_qual_points)
{
    // rather essential- copy the shared_ptr to the underlying mesh type!
    assert(mesh_lib_id <= mesh_library.size());
	mesh_ptr = boost::shared_ptr<ListedMesh>(mesh_library[mesh_lib_id]);
    
    init();
    set_dielectrics(protein_dielectric, solvent_dielectric);
}

void MeshInstance::init()
{
    const Mesh& mesh = *mesh_ptr;

    // set BasicNodePatches
    patches.clear();
    const std::vector<BasicNodePatch>& library_node_patches
		= mesh.get_node_patches();
    for (std::vector<BasicNodePatch>::const_iterator
			it=library_node_patches.cbegin(), end=library_node_patches.cend();
		 it != end; ++it)
    {
        // copy the node patches from the reference mesh
        // create the NodePatch but immediately store it as a pointer to the
        // base class;  the NodePatch constructor does change_coordinate_frame.
        // It is important that order matches library mesh.
        boost::shared_ptr<BasicNodePatch> ptr(new NodePatch(*it, *this)); 
        patches.push_back(ptr);
    }

    // set Charges
    charges.clear();
    const std::vector<Charge>& library_charges = mesh.get_charges(true);
    for (std::vector<Charge>::const_iterator
			it=library_charges.cbegin(), end=library_charges.cend();
		 it != end; ++it)
    {
		boost::shared_ptr<Charge> ptr(
			new Charge(*it, mesh.get_centre(), rotation, xyz_offset)
		);
		// Important: the same underlying charge is referenced twice!
		allCharges.push_back(ptr);	// Copy
		if ((*it).charge != 0) charges.push_back(std::move(ptr));	// Move
    }
    
    // set radius
    radius = mesh.get_radius();
}

static bool pt_triangle(const Vector& pt, const Vector* tv, int zdir,
						meshing::PtSet& hitVertices) {
	return Meshing<MeshInstance>::pt_triangle(pt, tv, zdir, hitVertices);
}

static bool pt_triangle(const Vector& pt, const Vector* tv, const Vector& dir,
						meshing::PtSet& hitVertices, Vector *pz) {
	return Meshing<MeshInstance>::pt_triangle(pt, tv, dir, hitVertices, pz);
}

bool MeshInstance::pt_is_internal(const Vector& pt,
		meshing::PtSet& hitVertices, const Vector* dir, Vector* nearest) const {
	return Meshing<MeshInstance>::pt_is_internal(*this, pt, hitVertices,
												 dir, nearest);
}

double MeshInstance::get_potential_at_internal_pt(const Vector& pt) const
{
    assert(this->pt_is_internal(pt) == true);

    // Assuming the point is internal to the mesh, the potential
    // is the sum of the coulomb terms from the charges and the
    // induced charge/dipole contribs from the surface.
    double pot=0.0;

    if (charge_fmm_tree.get() == nullptr) {
        charge_fmm_tree.reset
			(new fmm::FMM_Octree_6FIG_ACCURACY(10000, xyz_offset, radius * 4.));
        for (std::vector<boost::shared_ptr<Charge> >::iterator
				charge_it=charges.begin(), charge_end=charges.end();
             charge_it != charge_end; ++charge_it)
        {
            charge_fmm_tree->add(*charge_it);
        }
		// solve for kappa=0, i.e. non-screened Coulomb potential
        charge_fmm_tree->solve(0.0);
    }
    pot = charge_fmm_tree->calculate_potential(pt);
    pot *= ONE_OVER_4PI / Dprotein;

    // now work out the mesh integral terms
    // assume constant node patch treatment
    double accum=0.0;
    for (std::vector< boost::shared_ptr<BasicNodePatch> >::const_iterator
			np_it=patches.begin(), np_end=patches.end(); 
         np_it != np_end; ++np_it)
    {
        const BasicNodePatch& np = **np_it;
        const Vector& n = np.get_normal();

        double Gpt_val = Gpt(pt, np);
        double dGpt_val = dGpt_dn(pt, np, n);
        accum += np.get_bezier_area()
					* (-Gpt_val*np.h*np.dielectric_ratio + dGpt_val*np.f);
    }

    return pot+accum;
}

double MeshInstance::get_potential_at_external_pt(
	const Vector& pt,
	double kappa) const
{
    assert(this->pt_is_internal(pt) == false);

    // now work out the mesh integral terms
    // assume constant node patch treatment
    double accum=0.0;
    for (std::vector< boost::shared_ptr<BasicNodePatch> >::const_iterator
			np_it=patches.cbegin(), np_end=patches.cend(); 
         np_it != np_end; ++np_it)
    {
        const BasicNodePatch& np = **np_it;
        const Vector& n = np.get_normal();

        double upt_val = upt(pt, np, kappa);
        double dupt_val = dupt_dn(pt, np, n, kappa);
        accum += np.get_bezier_area() * (dupt_val*np.f - upt_val*np.h);
    }
    
    return accum;
}

// Rotation and translation is applied to a mesh instance not a reference mesh
MeshInstance&
	MeshInstance::move(const Vector& translate, const Quaternion& rotate)
{
std::cout << "MeshInstance::move from " << xyz_offset << " to " << translate << std::endl;
	for (auto& spbnp: patches) {
		// Change node patch coordinates
		(static_cast<NodePatch&>(*spbnp)).change_coordinate_frame
										(xyz_offset, rotate, translate);
		// Reset the potential, its derivative and non-electrostatic energies
		spbnp->f = 0.0;
		spbnp->h = 0.0;
		spbnp->he = 0.0;
		spbnp->lj = 0.0;
    }

    // set Charges:  relies on allCharges and charges pointing to same Charges
	for (auto& spch: allCharges) {
		spch->Vector::change_coordinate_frame(xyz_offset, rotate, translate);
    }
    
	// Update object copy of position and rotation - before resetting tree!!
    xyz_offset = translate;
    rotation *= rotate;

	return *this;
}

// Local values and private methods to calculate non-electrostatic energies
// Lambda to test if range is too great to be worth fiddling with
static auto toofar = [](Vector& v, double lim) noexcept -> bool {
	for (int i = 0; i < 3; i++) if (abs(v(i)) > lim) return true;
	return false;
};

static constexpr double EINT{10000};	// Internal point cost: max energy/area
static constexpr double P5olim{8.55};		// HE P5 Outer limit
//#define P5POLY
#ifdef P5POLY
static constexpr double P5ilim{3.64356};	// HE P5 Inner limit
#else  // P5POLY
static constexpr double P5mlim{6.8};		// HE P5 Middle limit
static constexpr double P5ilim{4.3};		// HE P5 Inner limit
#endif  // P5POLY
static constexpr double LJolim{2.5};		// LJ Outer limit
static constexpr double LJilim{2.3};		// LJ Inner limit
static constexpr double ELIM{std::max<double>(P5olim, 3.4*LJolim)};	//sigma<3.4 

// Useful debug tools
static std::mutex dbgmx;
#define COUT { std::lock_guard<std::mutex> dbglk(dbgmx); std::cout
#define ENDL std::endl; }

// Returns false if other not in range of mine under P5 (hydrophobic effect)
static bool P5range(const BasicNodePatch& mine, const BasicNodePatch& other) {
//COUT << "MeshInstance::P5 range in " << ENDL

	// No point doing anything if the hydrophobicity isn't right
	if (mine.hydrophobicity >= 0) return false;
	// positive values should have an effect through longer-range smoothing...
	
	// Cheap range checks
	Vector v = mine.get_node() - other.get_node();
	if (toofar(v, P5olim)) return false;

	// Calculate P5
	// Check range and approach
	double r = v.length();
	if (r >= P5olim) return false;		// too far away
	if (r > 0 && v.dot(mine.get_normal()) >= 0) return false; // opposite faces
	double n = mine.get_normal().dot(other.get_normal());
	if (n >= 0) return false;			// not enclosing water

//COUT << "MeshInstance::P5 range out " << ENDL
	// Otherwise, other is in range of mine for P5
	return true;
}

// Returns the P5 (hydrophobic "energy") for distance r and given hydrophobicity
static double P5(double r, double hydro) {
#ifdef P5POLY
	constexpr double l[] = {-0.1381608, 0.0395076 }; // Linear coefficients
	constexpr double q[] = {-3.256875, 2.82834, -0.957078, 0.158458,
							-0.0128423, 0.0004076924}; // Quintic coeffs
#else  // P5POLY
	constexpr double l[] = {-0.4144824, 0.1185228 }; // Linear coefficients
	constexpr double qm[] = {0.0245255, -4.0755e-3};  // Quadratic-middle coeffs
	constexpr double qo[] = {0.0261969, -3.0642e-3};  // Quadratic-outer coeffs
#endif  // P5POLY

	// Return value
	double e = 0.0;

//COUT << "MeshInstance::P5 in " << r << " " << hydro << ENDL
	// Apply the linear (inner) or quintic (outer) potentials
	if (r <= P5ilim) {
		//if (r < 0.0) r = 0.0;  // overlaps reduce volume but treat as contact
		e = l[1]*r + l[0];
	}
#ifdef P5POLY
	else // (P5ilim < r < P5olim)
		e = ((((q[5]*r + q[4])*r + q[3])*r + q[2])*r + q[1])*r + q[0];
#else  // P5POLY
	else if (r <= P5mlim)
		e = (qm[1]*r + qm[0])*r;
	else // P5mlim < r < P5olim
		e = (qo[1]*r + qo[0])*r;
#endif  // P5POLY

	// Adjust by strength of hydrophobicity - may have been smoothed here
	if (hydro != -0.5) e *= (-2)*hydro;
//COUT << "MeshInstance::P5 out " << e << ENDL
	return e;
}

static bool LJrange(const Charge& mine, const Charge& other) {
	// No point doing anything if parameters are set 0
//COUT << "MeshInstance::LJ range in " << ENDL
	if (other.epsilon == 0.0 || other.sigma == 0.0) return false;

	// Cheap range checks
	Vector v = mine - other;

	if (toofar(v, LJolim*other.sigma)) return false;

	// In this case it is more efficient to check the outer limit exactly
	// during the LJ call, so this is deferred

//COUT << "MeshInstance::LJ range out " << ENDL
	// Otherwise, other is in range of mine for LJ
	return true;
}

static double LJ(double r, double sigma, double epsilon) {
	//constexpr double m[] = {3136.5686, -68.069, -0.0833111261, 0.746882273};
	//constexpr double k = 0.0163169237;
	constexpr double m[] = {1568.2843, -34.0345, -0.041655563, 0.373441136};
	constexpr double k = 0.0081584618;

//COUT << "MeshInstance::LJ in " << r << " " << sigma << " " << epsilon << ENDL
	// Check limits
	if (r >= LJolim*sigma) return 0.0;
	if (r == 0.0) return EINT;

	// Calculate LJ energy
	double sr2 = sigma/r;
	sr2 *= sr2;
	double sr6 = sr2*sr2*sr2;
	double sr12 = sr6*sr6;

	double e; // return value
	if (r > LJilim*sigma)	// (LJilim < r < LJolim)
		e = (m[0]*sr12 + m[1]*sr6 + m[2]/sr2 + m[3]) * epsilon;
	else // (0.0 < r <= LJilim sigma)
		e = (2 * (sr12 - sr6) + k) * epsilon;
	if (e > EINT) e = EINT;

//COUT << "MeshInstance::LJ out " << e << ENDL
	return e;
}

unsigned int MeshInstance::resize_pot(unsigned int n) {
	if (n > hepot[psel].size()) {
		for (auto& pot: {hepot, ljpot}) {
			for (int s = 0; s < 2; s++) pot[s].resize(n); // init 0.0 extras
		}
	}
	return n;
}

// Local types and class to handle threading of calculate_boundary_energy
// Kahan summation to maintain totals for mesh instance
struct kahanValue {
	double sum, carry;
	kahanValue() { sum = 0.0; carry = 0.0; };
};
template<unsigned int N> struct kahanValueList { kahanValue value[N]; };

static auto kahan = [](double& t, double& c, double el) mutable -> double {
	el -= c;				// Adjusted element
	double s = t + el;		// Imprecise sum
	c = (s - t) - el;		// Carried error
	t = s;					// Store the result
	return t;
};

static auto kahanSum = [](kahanValue& kv, double el) mutable -> double {
	return kahan(kv.sum, kv.carry, el);
};

// Utility to generate a random list of integers [0, n)
static std::vector<unsigned int>shuffle(unsigned int n) {
	std::vector<unsigned int> rv(n);
	std::iota(rv.begin(), rv.end(), 0);
	std::shuffle(rv.begin(), rv.end(), std::mt19937{std::random_device{}()});
	return rv;
}

// thread class for calculate_boundary_energy, but actually a little more
// generic  e.g. supports work queuing although a single constructor list
// would have done.
struct CBE_thread {
	// Primary (only) constructor:
	// result contains 0: energy, 1: area
	// trklen is the number of "omi" patches
	CBE_thread(kahanValueList<3>& result, unsigned int trklen) :
		gen(rd()),
		uniform(0.0, 1.0)	// Uniform [0,1)
	{
		trk.resize(trklen);
		retval = &result;
		waiting = false;	// See below, suggests a code issue
		th = std::thread(&CBE_thread::process, this);
	};
	CBE_thread(const CBE_thread&) = delete;	// Not copyable
	CBE_thread(CBE_thread&&) = default;		// But is movable
	~CBE_thread(void) { stop(); };			// Implied stop

	// Explicit stop:  notify all and wait for thread to finish its work
	void stop(void) {
		bool already_stopping = CBE_thread::stopping;
		CBE_thread::stopping = true;
		if (!already_stopping) CBE_thread::cv.notify_all();
		// This ought not to be necessary but it appears that
		// a) without the lock worker thread can be in between checking
		// condition and unlock+wait within cv.wait(), thereby missing both
		// the change to stopping and the notification.
		// b) without the waiting flag, there is a deadlock
		// Probably the code is not structured optimally...
		if (waiting) {
			std::lock_guard<std::mutex> lock(CBE_thread::lkq);
			cv.notify_all();
		}
		th.join();
	}
	// Immediate stop: a reset is required after this
	static void halt(void) { halting = true; };

	// Queue a work list
	using cbeFn = std::function< void(
		BasicNodePatch&, std::vector<unsigned int>&, meshing::PtSet&) >;
	using cbeTpl = std::tuple<cbeFn, PatchList::iterator, PatchList::iterator>;
	//static void push(const cbeTpl& work_item) {	// no std::forward<T> below
	static void push(cbeTpl&& work_item) {
		{	// new scope
			std::lock_guard<std::mutex> lock(CBE_thread::lkq);
			CBE_thread::q.push(std::make_shared<cbeTpl>(
									std::forward<cbeTpl>(work_item)));
		}
		CBE_thread::cv.notify_one();
	};
	// Reset queue and flags
	static void reset(void) {
		std::queue<std::shared_ptr<cbeTpl>> empty;
		CBE_thread::q.swap(empty);
		CBE_thread::stopping = false;
		CBE_thread::halting = false;
	}

private:
	// Class attributes to coodinate all threads
	static std::queue<std::shared_ptr<cbeTpl>> q;
	static std::mutex lkq;
	static std::condition_variable cv;
	static std::atomic<bool> stopping, halting;

	// These are needed to support randomisation
	std::random_device rd;		// obtain seed for random number engine
	std::mt19937_64 gen;		// 64-bit mersenne_twister_engine
	std::uniform_real_distribution<double> uniform;

	// Object attributes
	std::thread th;
	kahanValueList<3> *retval;
	std::atomic<bool> waiting;

	// Tracking objects
	// This could be generalised by templating
	meshing::PtSet hitVertices;
	std::vector<unsigned int> trk;

	// Main processing
	void process(void) {
		while (!CBE_thread::stopping || !CBE_thread::q.empty()) {
			cbeFn work;  // zero-initialised
			PatchList::iterator begin, end;
			{//new scope
				std::unique_lock<std::mutex> lk(CBE_thread::lkq);
				waiting = true;
				CBE_thread::cv.wait(lk, []{return CBE_thread::stopping ||
												!CBE_thread::q.empty();});
				waiting = false;
				if (CBE_thread::stopping && CBE_thread::q.empty()) break;
				if (!CBE_thread::q.empty()) {
					std::tie(work, begin, end) = *CBE_thread::q.front();
					CBE_thread::q.pop();
				}
			}
			//CBE_thread::cv.notify_one();	// Not needed! All waiters check q
			if (work) {
				for (auto&& it = begin; it != end; it++) {
					if (halting) break;
					BasicNodePatch& np = **it;

					// It is possible there are multiple mesh instances in
					// range of this one - select later ones through MC
					double che = np.he;
					static const double kT = beep_constants::Boltzmann_SI
											* beep_constants::Avogadro / 1000;

					// Now run the cbe method
					// This could be generalised by returning a tuple to tie
					work(np, trk, hitVertices);

					// Update totals - generalise by moving into cbe
					// For HE, use a simple MC to decide which to keep
					if (che == 0.0 || np.he <= che ||
						exp((che-np.he)/kT) > uniform(gen)) {
						// Accepting the new value - adjust current total
						kahanSum(retval->value[0], np.he-che);
					}
					// Rejecting the new one - reset the NP value
					else
						np.he = che;

					// LJ is just summed
					kahanSum(retval->value[1], np.lj);

					// Not perfect, but try to avoid counting area twice
					if (!np.hea && np.he != 0.0) {
						double a = np.get_bezier_area();
						kahanSum(retval->value[2], a);
						np.hea = true;
					}
				}
			}
		}
	};
};
std::queue<std::shared_ptr<CBE_thread::cbeTpl>> CBE_thread::q;
std::mutex CBE_thread::lkq;
std::condition_variable CBE_thread::cv;
std::atomic<bool> CBE_thread::stopping, CBE_thread::halting;

void MeshInstance::calculate_boundary_energy(
	BasicNodePatch& np,						// The patch to calculate for
	const MeshInstance& omi,				// The impinging MeshInstance
	std::vector<unsigned int>& track,		// To avoid repeat omi points
	meshing::PtSet& hitVertices)			// To avoid repeat vertex hits
{
	// Calculate the energies for np resulting from omi
	// Reset the energies
	np.he = 0.0;
	np.lj = 0.0;

	// Get the point and the list of omi triangles
    const Vector& pt = np.get_node();
	const auto& triangles = omi.mesh_ptr->get_triangle_ptrs();

	// PART 1:  Calculate provisional energy and determine if pt is inside omi
	// Set up pt_triangle inputs
	Vector tv[3];	// Triangle - copies because const pointers are a mess
	const int zdir = 2*(pt.z - omi.xyz_offset.z > 0.0)-1;	// z-dir to use +/-
	hitVertices.clear();	// See interface to pt_triangle
	int i;			// General index
	int cc = 0;		// crossings count for z-line

	unsigned int npatches = omi.patches.size();	// = number of omi vertices
	unsigned int start = 0;		// To track secondary encounters of omi vertices
	if (npatches > 0) start = track[0];	// All values should be the same

	// Locate this patch's charge, which will be needed for LJ potential
	const Charge *mych = nullptr;
	if (np.ch_idx < allCharges.size()) mych = allCharges[np.ch_idx].get();

	// Set up the unordered sets used to hold vertices, triangles and charges
	// These should be quite small, so can be set up for each call
	using uint = unsigned int;
	using idxSet = std::unordered_set<uint>;
	idxSet vertexSet;
	using idxMap = std::unordered_map<uint, bool>;
	idxMap chargeMap;
	using trio = std::array<uint,3>;
	struct hashFnTrio {
		size_t operator()(const trio& v) const{ return boost::hash_value(v); }
	};
	struct eqFnTrio  {
		bool operator()(const trio& lhs, const trio& rhs) const
		{ return (lhs == rhs); }
	};
	
	using idxTri = std::unordered_set<trio, hashFnTrio, eqFnTrio>;
	idxTri triangleSet;

	// Loop over all omi triangles checking potential ranges whilst
	// also seeking crossings in z-dir from pt
    for (const auto& spt: triangles) {
		const Triangle& t = *spt;
		std::array<uint,3> vi
			= {t.get_v1_idx(), t.get_v2_idx(), t.get_v3_idx()};
		assert(vi[0] < npatches && vi[1] < npatches && vi[2] < npatches);
		// This relies on vertices and patches having the same order...
		for (i = 0; i < 3; i++) {
			tv[i] = omi.patches[vi[i]]->get_node();

			// Calculate provisional boundary energy
			// Have we been here before?
			if (track[vi[i]] == start) {
				// Range checking for new vertex
				const BasicNodePatch& onp = *(omi.patches[vi[i]]);
				if (P5range(np, onp)) vertexSet.insert(vi[i]);

				// Locate the two charges, which contain LJ parameters
				// No charges means no LJ energy to be included
				if (mych != nullptr && onp.ch_idx < omi.allCharges.size()) {
					const Charge& otch = *(omi.allCharges[onp.ch_idx]);
					if (LJrange(*mych, otch))
						chargeMap.insert({onp.ch_idx, false});
				}
				track[vi[i]]++;
			}
			// add triangle to list if vertex is there
			if (vertexSet.count(vi[i]) > 0) triangleSet.insert(vi);
		}

		// Does z-line does cross the mesh at this triangle?
		cc += pt_triangle(pt, tv, zdir, hitVertices);

	}  // each triangle in mesh

	// pt is internal to omi if odd number of crossings
	bool internal = false;
	if ((cc%2) == 1) internal = true;

	// Hydrophobic effect energy - done first because has longer range
	// and can cut out if still too far.
	double r = -1.0;	// distance for P5, negative means HE is 0
	double n = -1.0;	// normal.other-normal, so must be negative
	// Cut out now if nothing was in range
	if (vertexSet.size() != 0) {
		// Not internal and not distant, so calculate the energy
		// The hydrophobic effect requires a distance, which is in the mean
		// direction to the omi patches in range from pt.
		unsigned ctr = 1;
		Vector dir(0,0,0);
		const Vector origin(0,0,0);
		for (const auto& idx: vertexSet) {
			Vector v = omi.patches[idx]->get_node() - pt;
			if (v != origin) dir += (v.normalised() - dir)/ctr++;
		}
		// Last resort - this should never happen
		if (dir == origin) dir = np.get_normal();

		// Find the minimum distance to an in-range triangle in this direction.
		// This will be on the near side of the omi.
		double rmin = ELIM;
		unsigned int vidx = 0;  // arbitrary choice if nothing beats ELIM
		Vector opt;
		for (const auto& aidx: triangleSet) {
			for (i = 0; i < 3; i++) tv[i] = omi.patches[aidx[i]]->get_node();
			Vector ptmp;
			if (pt_triangle(pt, tv, dir, hitVertices, &ptmp)) {
				for (unsigned int i = 0; i < 3; i++) {
					double rtmp = (tv[i]-ptmp).length();
					if (rtmp < rmin) {
						vidx = aidx[i];
						opt = ptmp;
						rmin = rtmp;
					}
				}
			}
		}
		if (rmin >= ELIM) return;	// Turned out nothing in range - how?

		// Now calculate the energy per unit area for this local direction
		// Convert energy per unit area to energy and adjust by node angle
		const BasicNodePatch& onp = *(omi.patches[vidx]);
		r = (pt-opt).length();
		n = np.get_normal().dot(onp.get_normal());  // will be negative
	}
	if (r >= 0) {
		double a = np.get_bezier_area();	// seems planar if force_planar
		np.he = P5(r*(internal? -1: 1), np.hydrophobicity) * a*n*(-1);
		if (np.he > EINT) np.he = EINT;		// should not happen for P5
	}

	// For Lennard-Jones, assume additivity of dispersion and repulsion energies
	// Locate the two charges, which contain most of the required information
	if (mych != nullptr) {
		for (const auto& pr: chargeMap) {
			const Charge& otch = *(omi.allCharges[pr.first]);
			if (internal && !pr.second) {
				chargeMap[pr.first] = true;	// Show we have been here before
				// This is the first time this charge has been tested for this
				// internal node - this is an expensive test...
				hitVertices.clear();
				if (omi.pt_is_internal(*mych, hitVertices)) {
					// Bad news - EINT and halt all other processing
					//CBE_thread::halt();	// Continue to probe HE energy
					//np.he = 0.0;	// almost certainly unnecessary
					np.lj = EINT;
					std::cout << "LJ collision detected between mesh instances "
							  << instance_id << " and "
							  << omi.instance_id << std::endl;
					return;
				}
				// Good news - processing can continue as normal...
			}
			double r = (*mych - otch).length();
			np.lj += LJ(r, otch.sigma, otch.epsilon) /
					mesh_ptr->get_npcount(np.ch_idx);
		}
		if (np.lj > EINT) np.lj = EINT;
	}
}

void MeshInstance::update_energy(MeshInstanceList& mis) {
	static unsigned int num_threads = (std::thread::hardware_concurrency() > 0 ?
		std::thread::hardware_concurrency() : 1);
	std::cout << num_threads << " threads" << std::endl;

	// Hydrophobic effect is only valid relative to the moving MI, and
	// area of contact is only the new contact area for the move
	hea = 0.0;
	for (auto& nit: patches) { nit->he = 0.0; nit->hea = false; }

	// Update the energy resulting from all other MIs
	for (auto& spmi: mis) {	// no return or break;  continue is ok
		MeshInstance& omi = *spmi;
		if (omi.instance_id == instance_id) continue;
		resize_pot(omi.instance_id+1);
		omi.resize_pot(instance_id+1);

		// Save current values as now previous, reset current
		MeshInstance* m[2] = {this, &omi};
		for (int i = 0; i < 2; i++) {
			unsigned int oid = m[1-i]->instance_id;
			for (auto& pot: {m[i]->hepot, m[i]->ljpot}) {
				pot[1-m[i]->psel][oid] = pot[m[i]->psel][oid];
				pot[m[i]->psel][oid] = 0.0;
			}
		}
		// And reset the he data for this MI
		omi.hea = 0.0;
		for (auto& onit: omi.patches) { onit->he = 0.0; onit->hea = false; }

		// Fast check on range: too far apart and no possible interaction
		if ((xyz_offset - omi.xyz_offset).length() > radius + omi.radius + ELIM)
			continue;

		// Update energy from impact of other mesh instances
		MeshInstance* mi[2] = {this, &omi};
		for (int i = 0; i < 2; i++) {	// no return or break;  continue is ok
			auto& nps = mi[i]->patches;
			auto nit = nps.begin();
			unsigned int npatches = nps.size();
			unsigned int onpatches = mi[1-i]->patches.size();
			std::vector<kahanValueList<3>> r(num_threads);

			// Start threads in a new scope to limit their existence
			if (num_threads > 1) {
				// make sure boost isn't messing up placeholders!!!
				using std::placeholders::_1; using std::placeholders::_2;
				using std::placeholders::_3;

				CBE_thread::reset();	// clear queue and unset flags
				std::vector<std::unique_ptr<CBE_thread>> thl;
				for (int j = num_threads; j > 0; j--) {
					thl.push_back(
						std::make_unique<CBE_thread>(r[j-1],onpatches) );
					auto f = std::bind(&MeshInstance::calculate_boundary_energy,
										mi[i], _1, *(mi[1-i]), _2, _3);
					unsigned int block = npatches/j;
					auto nite = (npatches > block ? nit + block : nps.end());
					CBE_thread::push(std::make_tuple(f, nit, nite));

					npatches -= block;
					if (npatches <= 0) break;
					nit += block;
				}
			}	// Main thread sits here until worker threads complete
			else {
				// threading not in vogue - just run in main thread
				std::vector<unsigned int>trk(onpatches, 0);
				static const double kT = beep_constants::Boltzmann_SI
											* beep_constants::Avogadro / 1000;
				// These are needed to support randomisation
				std::random_device rd;
				std::mt19937_64 gen(rd());
				std::uniform_real_distribution<double> uniform(0.0, 1.0);

				for (auto& nit: nps) {
					BasicNodePatch& np = *nit;
					// It is possible there are multiple mesh instances in
					// range of this one - select later ones through MC
					double che = np.he;
					calculate_boundary_energy(np, *(mi[1-i]), trk);

					// Update totals
					// For HE, use a simple MC to decide which to keep
					if (che == 0.0 || np.he < che ||
						exp((che-np.he)/kT) > uniform(gen)) {
						// Accepting the new value - adjust current total
						kahanSum(r[0].value[0], np.he-che);
					}
					// Rejecting the new one - reset the NP value
					else
						np.he = che;

					// LJ is just summed
					kahanSum(r[0].value[1], np.lj);

					// Avoid counting area twice
					if (!np.hea && np.he != 0.0) {
						double a = np.get_bezier_area();
						kahanSum(r[0].value[2], a);
						np.hea = true;
					}
				}
			}

			// Final sums and limit checks (threads completed above)
			double ce[2] = {0.0, 0.0}, ca = 0.0;
			for (auto& kv: r) {
				int val = 0;
				for (auto mipot:
						{&mi[i]->hepot[mi[i]->psel][mi[1-i]->instance_id],
						 &mi[i]->ljpot[mi[i]->psel][mi[1-i]->instance_id]})
				{
					ce[val] += kv.value[val].carry;
					kahan(*mipot, ce[val], kv.value[val].sum);
					val++;
					if (*mipot > EINT) *mipot = EINT;
				}
				ca += kv.value[val].carry;
				kahan(mi[i]->hea, ca, kv.value[val].sum);
			}
		}	// each way

		// Report the changes in total and hydrophobic areas
		std::cout << "Change in hydrophobic area of contact instance "
				<< omi.instance_id << " is " << omi.hea << std::endl;
	}	// each other MI
	std::cout << "Change in hydrophobic area of contact instance "
			<< instance_id << " is " << hea << std::endl;
}

void MeshInstance::revert_energy(MeshInstanceList& mis) {
	for (auto& spmi: mis) {
		if (spmi->instance_id == instance_id) continue;
		for (auto& pot: {spmi->hepot, spmi->ljpot})
			pot[spmi->psel][instance_id] = pot[1-spmi->psel][instance_id];
	}
	psel = 1-psel;
}

double MeshInstance::calculate_energy(
	bool electrostatics,
	double kappa,
	double fvals[],
	double hvals[]) const
{
    // Dsolvent and Dprotein should have been set via set_dielectrics()
	double E = 0.0;
	if (electrostatics)
		E = mesh_ptr->calculate_energy(kappa, Dprotein, Dsolvent,
	                                      fvals, hvals);
    std::cout << "Energy for mesh " << instance_id << " (lib_id="
	          << mesh_ptr->get_id() << ")" << std::setprecision(10)
			  << " electrostatic=" << E;
	double c = 0.0, s;
	const std::vector<double>* pot[] = {hepot, ljpot};
	const char* name[] = {" HE=", " LJ="};
	for (int p = 0; p < 2; p++) {
		double F = E;
		for (const auto& e: pot[p][psel]) {
			// Kahan sum of contributions
			s = E + (e - c);
			c += (s - E) - e;
			E = s;
		}
		std::cout << name[p] << E-F;
	}
	std::cout << " total=" << E << std::endl;
    return E;
}

Vector MeshInstance::calculate_force(
	double kappa,
	double fvals[],
	double hvals[]) const
{
    // Dsolvent and Dprotein should have been set via set_dielectrics()
    KahanVector MST_ext,MST_int,dbf,ionic;
    mesh_ptr->calculate_surface_integral_forces(kappa, Dprotein, Dsolvent,
								 fvals, hvals, MST_ext, MST_int, dbf, ionic);
    Vector force(*MST_ext);
    force.apply_rotation(rotation);
    return force;
}

void MeshInstance::calculate_forces(
	double kappa,
	double fvals[],
	double hvals[],
	KahanVector& qE,
	KahanVector& MST_ext,
	KahanVector& MST_int,
	KahanVector& dbf,
	KahanVector& ionic) const
{
    // Dsolvent and Dprotein should have been via the set_dielectrics()
    qE = mesh_ptr->calculate_qe_force(Dprotein, Dsolvent, fvals, hvals);
    mesh_ptr->calculate_surface_integral_forces(kappa, Dprotein, Dsolvent,
								 fvals, hvals, MST_ext, MST_int, dbf, ionic);
    
    (qE).apply_rotation(rotation);
    (MST_int).apply_rotation(rotation);
    (MST_ext).apply_rotation(rotation);
    (ionic).apply_rotation(rotation);
    (dbf).apply_rotation(rotation);
}

	static constexpr unsigned int kin_nvals = 7;

void MeshInstance::kinemage_fh_vals(
	double fscale,
	double hscale,
	int num_colours,
	std::ostringstream& buf_f,
	std::ostringstream& buf_h) const
{
	double s[kin_nvals] = {0.0, 0.0, 0.0, 0.0, 0.0, fscale, hscale};
	std::ostringstream buf_hy, buf_s, buf_e, buf_he, buf_lj;
	kinemage_vals(num_colours, buf_hy, buf_s, buf_e, buf_he, buf_lj,
				  buf_f, buf_h, s);
}

void MeshInstance::kinemage_vals(
	int num_colours,
	std::ostringstream& buf_hy,
	std::ostringstream& buf_s,
	std::ostringstream& buf_e,
	std::ostringstream& buf_he,
	std::ostringstream& buf_lj,
	std::ostringstream& buf_f,
	std::ostringstream& buf_h) const
{
	double s[kin_nvals]{};
	kinemage_vals(num_colours, buf_hy, buf_s, buf_e, buf_he, buf_lj,
				  buf_f, buf_h, s);
}

template<unsigned int S>
void MeshInstance::kinemage_vals(
	int num_colours,
	std::ostringstream& buf_hy,
	std::ostringstream& buf_s,
	std::ostringstream& buf_e,
	std::ostringstream& buf_he,
	std::ostringstream& buf_lj,
	std::ostringstream& buf_f,
	std::ostringstream& buf_h,
	double (&scale)[S]) const
{
	if (S != kin_nvals) {
		std::cerr << "MeshInstance::kinemage_vals invalid scales"
				  << std::endl;
		throw std::exception();
	}

	// Obtain scales
	unsigned int npatches = 0;
	bool getscale[kin_nvals];
	for (int i = 0; i < kin_nvals; i++) {
		getscale[i] = (scale[i] == 0.0);
		if (getscale[i]) npatches = patches.size();
	}
    for (size_t np_ctr = 0; np_ctr < npatches; ++np_ctr) {
        const NodePatch& np = dynamic_cast<NodePatch&>(*(patches[np_ctr]));
		double sig = 0.0, eps = 0.0;
		if (np.ch_idx < allCharges.size()) {
			Charge& ch = *(allCharges[np.ch_idx]);
			sig = ch.sigma;
			eps = ch.epsilon;
		}
		double value[] = {np.hydrophobicity, sig, eps,
						  np.he, np.lj, np.f, np.h};
		for (int i = 0; i < kin_nvals; i++)
			if (getscale[i])
				scale[i] = std::max<double>(scale[i], fabs(value[i]));
	}
	for (int i = 0; i < kin_nvals; i++)
{
		if (scale[i] == 0.0) scale[i] = 1.0;
std::cout << __PRETTY_FUNCTION__ << " scale " << i << " is " << scale[i] << std::endl;
}

	// Output colours
    for (size_t np_ctr = 0; np_ctr < patches.size(); ++np_ctr) {
        const NodePatch& np = dynamic_cast<NodePatch&>(*(patches[np_ctr]));
        
		double sig = 0.0, eps = 0.0;
		if (np.ch_idx < allCharges.size()) {
			Charge& ch = *(allCharges[np.ch_idx]);
			sig = ch.sigma;
			eps = ch.epsilon;
		}
		double value[] = {np.hydrophobicity, sig, eps,
						  np.he, np.lj, np.f, np.h};
		std::ostringstream *buf[] = {&buf_hy, &buf_s, &buf_e,
									 &buf_he, &buf_lj, &buf_f, &buf_h};
        std::string name[kin_nvals];

        // figure out the colour name corresponding to the values
        for (int i = 0; i < kin_nvals; i++) {
            std::stringstream s;
            int idx = static_cast<int>(round(
					static_cast<double>(num_colours)*value[i]/scale[i]));
            if (idx < 1){
                idx = abs(idx);
                int col = (idx <= num_colours ? idx : num_colours+1);
                s << "red_" << col;
            } else {
                int col = (idx <= num_colours ? idx : num_colours+1);
                s << "blue_" << col;
if (col > num_colours) {
std::cout << __PRETTY_FUNCTION__ << " excess " << i << " " << value[i] << " " << scale[i] << std::endl;
}
            }
            name[i] = s.str();
        }

        Vector local_centre = mesh_ptr->get_centre();
        std::vector<PointNormal> edge_points;
        np.get_edge_points(edge_points);
        for (auto it = edge_points.cbegin(); it != edge_points.cend(); ++it) {
            auto next_it = it+1;
            if (next_it == edge_points.cend()) next_it = edge_points.begin();
            
            Vector here = it->pt();
            here.change_coordinate_frame(local_centre, rotation, xyz_offset);
            Vector next = next_it->pt();
            next.change_coordinate_frame(local_centre, rotation, xyz_offset);

            // kinemage a quadrilateral with the values
            for (int i = 0; i < kin_nvals; i++) {
				*buf[i] << "X " << name[i] << " "
						<< np.x << " " << np.y << " " << np.z << " "
						<< name[i] << " "
						<< here.x << " " << here.y << " " << here.z << " "
						<< name[i] << " "
						<< next.x << " " << next.y << " " << next.z << "\n";
			}
        }
    }
}

void MeshInstance::set_unique_patch_id(unsigned int patch_ctr)
{
	for (PatchList::iterator nit=patches.begin(), nend=patches.end();
		 nit != nend; ++nit)
	{
		(**nit).set_idx(patch_ctr++);
	}
}

bool MeshInstance::init_fh_vals(
	double *x,
	unsigned int& xctr,
	unsigned int offset,
	bool preconditioned) const
{
	const Mesh& mesh = *mesh_ptr;
	const std::vector<BasicNodePatch>& patches = mesh.get_node_patches();
	for (unsigned int ctr=0; ctr < patches.size(); ++ctr) {
		double f_preset = patches[ctr].f;
		double h_preset = patches[ctr].h;
		preconditioned = preconditioned || f_preset != 0.0 || h_preset != 0.0;
		x[xctr] = f_preset;
		x[xctr+offset] = h_preset;
		xctr++;
	}
	return preconditioned;
}

//NB where do f_lhs etc really belong?
unsigned int MeshInstance::reset_fh_vals(
	const boost::shared_array<double>& f_lhs, 
	const boost::shared_array<double>& h_lhs,
	unsigned int start_ctr)
{
	unsigned int inc_ctr = start_ctr;
	Mesh& mesh = *mesh_ptr;
	std::vector<BasicNodePatch>& patches = mesh.get_node_patches();
	for (unsigned int ctr=0; ctr < patches.size(); ++ctr) {
		BasicNodePatch& np = patches[ctr];
		np.f = f_lhs[inc_ctr];
		np.h = h_lhs[inc_ctr];
		++inc_ctr;
	}
	return inc_ctr;
}


// MeshInstanceList

// Add a new mesh instance to the list
boost::shared_ptr<MeshInstance>
MeshInstanceList::add(
	unsigned int mesh_lib_id,
	unsigned int mesh_instance_id,
	const Vector& offset,
	const Quaternion& rotation,
	double Dprotein, double Dsolvent,
	unsigned int num_quad_points,
	unsigned int num_qual_points,
	bool _silent)
{
	boost::shared_ptr<MeshInstance> ptr(
		new MeshInstance(mesh_lib_id, mesh_instance_id, library,
							offset, rotation, Dprotein, Dsolvent,
							num_quad_points, num_qual_points, _silent)
	);

//NB Comments inherited from BEEP::insert_mesh_instance
// minimum separation 5A between node patches
//     const double minimum_separation = 5.;
//     
//     // loop over the node patch points and ensure that they are not within   an angstrom
//     // of any other node patch already in the system
//     for (std::vector< boost::shared_ptr<MeshInstance> >::const_iterator      minst_it = begin(), minst_end=end(); minst_it != minst_end; ++    minst_it)
//     {
//         const MeshInstance& minst = **minst_it;
//         bool further_check_required = false;
//         if ((minst.xyz_offset - ptr->xyz_offset).length() <    (ptr->radius + minst.radius + minimum_separation)) {
//             further_check_required = true;
//         }
// 
//         if (further_check_required)
//         {
// 
//             throw std::exception();
//             
// /*            for (PatchList::const_iterator
//                 
//                 
//             for (PatchList::const_iterator np_it = minst.patches.begin(),    np_end=minst.patches.end(); np_it != np_end; np_it != np_end)
//             {
//                 const BasicNodePatch& np = **np_it;
//                 np
//             }*/
//             
//         }
//         
//     }

	// set unique id on each node patch
	unsigned int patch_ctr = get_total_patches();
	ptr->set_unique_patch_id(patch_ctr);

	push_back(ptr);
	return ptr;
}

// Move an instance in the list
// This exists because the simpler MeshInstance::move() is not working:
// the energy calculation for the moved instance is too high and it is not
// apparent why this should be.  This approach mimics clear and insert
// which appear to work individually.
// Note: cannot just copy construct the new MeshInstance for the same reason.
boost::shared_ptr<MeshInstance>
MeshInstanceList::move(
	unsigned int mesh_instance_id,
	const Vector& offset,
	const Quaternion& rotation,
	double Dprotein, double Dsolvent,	//TODO NB Passed because not accessible!
	bool _silent)						// Usage?? Can they be removed?
{
	// Validate the id
    if (mesh_instance_id >= size()) {
        std::cerr << "Bad value for mesh instance_id: " << mesh_instance_id
			<< " (" << size() << " mesh instances defined)" << std::endl;
        throw std::exception();
    }
{
static bool warn = true;
if (warn)
std::cerr << "MeshInstanceList::move NEEDS REVIEW - will not work with hydro"
			<< std::endl;
warn = false;
}
	
	// Obtain require values from existing MeshInstance
	MeshInstance& m = *(this->at(mesh_instance_id));
	unsigned int mesh_lib_id
		= static_cast<const ListedMesh&>(m.get_ref_mesh()).get_id();
	unsigned int num_quad_points = m.get_quad_points_per_triangle();
	unsigned int num_qual_points = m.get_qual_points_per_triangle();
	unsigned int patch_ctr = (**(m.get_node_patches().begin())).get_idx();
	
	// Create a new instance
	boost::shared_ptr<MeshInstance> ptr(
		new MeshInstance(mesh_lib_id, mesh_instance_id, library,
							offset, rotation, Dprotein, Dsolvent,
							num_quad_points, num_qual_points, _silent)
	);

	// set unique id on each node patch
	ptr->set_unique_patch_id(patch_ctr);

	// Insert the new instance at the required point
	erase(begin()+mesh_instance_id);
	insert(begin()+mesh_instance_id, ptr);
	return ptr;
}


// INIT INITIAL VALUES to what was set in input fh files
bool MeshInstanceList::init_library_fh_vals(double *x, unsigned int offset)
{
    bool preconditioned = false;
    unsigned int xctr=0;
    for (auto& spmi: *this) {
		preconditioned = preconditioned
		               || spmi->init_fh_vals(x, xctr, offset, preconditioned);
    }
	return preconditioned;
}

// Reset the f/h node patch values of the mesh definitions
void MeshInstanceList::reset_library_fh_vals(
	const boost::shared_array<double>& f_lhs,
	const boost::shared_array<double>& h_lhs)
{
    unsigned int total_ctr=0;
    for (auto& spmi: *this)
		total_ctr += spmi->reset_fh_vals(f_lhs, h_lhs, total_ctr);
}

size_t MeshInstanceList::get_total_patches() const
{
    unsigned int total_np=0;
    for (const auto& spmi: *this)
        total_np += spmi->get_num_node_patches();
    return total_np;
}

