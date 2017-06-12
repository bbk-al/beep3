/*      Author: david fallaize   Created on: 21 Jul 2010 */
/*      Modified: adam light    on: 9 Mar 2012  */

/*! \file charge.h
 * \brief This module declares the Charge class.
 *	
 *	This module contains the following public classes:
 *	- Charge -- representing a partial charge within a mesh as a model of a
 *	molecule.
 */

#ifndef CHARGE_H_
#define CHARGE_H_

#include "../common/math_vector.h"
//#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <limits>
#include <string>
#include <vector>
#ifdef __CHARMC__
#include <pup_stl.h>
#endif

namespace fs = boost::filesystem;	// Easier to swap to std::filesystem in '17

class Charge : public Vector
{
public:
	// static methods
	inline static double read_charges_from_file(
		const fs::path& xyzqr_filename,
		std::vector<Charge>& charges,
		bool skip_zeroes=true);

    Charge() = default;							//! Default constructor
	Charge(const Charge&) = default;			//! Copy constructor
	Charge(Charge&&) = default;					//! Move constructor
	Charge& operator=(const Charge&) = default;	//! Copy assignment
	Charge& operator=(Charge&&) = default;		//! Move assignment
    ~Charge() = default;						//! Destructor

	//! Constructors from vector and charge
    Charge(const Vector& v, double ch) : Vector(v), charge(ch), radius(0) {}
    Charge(const Vector& v, double ch, double r)
	: Vector(v), charge(ch), radius(r) {}

	//! Variant copy constructor with translation and rotation applied
    Charge(const Charge& other,
           const Vector& centre_of_rotation,
           const Quaternion& rot,
           const Vector& xyz_offset) :
           Vector(other),
           charge(other.charge),
           radius(other.radius)
    {
        Vector::change_coordinate_frame(centre_of_rotation, rot, xyz_offset);
    }

	// get_ methods
    Vector& position() { return static_cast<Vector&>(*this); }
    const Vector& position() const { return static_cast<const Vector&>(*this); }
    //const Vector& py_get_position() const {
	//	return static_cast<const Vector&>(*this);
	//}
    double get_charge() const { return charge; }
    double& get_charge() { return charge; }
    double get_radius() const { return radius; }
    double& get_radius() { return radius; }

	//NB Argument is ignored??
    double get_charge(int) const { return charge; }
    double& get_charge(int) { return charge; }

#ifdef __CHARMC__

    /// PUP Routine ///
    void pup(PUP::er &p) {

        Vector::pup(p);
        p | charge;
        p | radius;

    }
#endif

// These are public for access from Python only.  Use get_* elsewhere.

    double charge;
    double radius;
};

typedef Charge CharmChargeHolder;

//NB LAZY:  this is here only to avoid creating a charge.cpp!
inline double Charge::read_charges_from_file(
	const fs::path& xyzqr_filename,
	std::vector<Charge>& charges,
	bool skip_zeroes)
{
	// open the file
	fs::ifstream xyzqr_file;
	xyzqr_file.open(xyzqr_filename);
	assert(xyzqr_file.good());

	double net_charge = 0.0;

	// read in xyzr data, and create charges
	double x,y,z,charge,rad;
	while(xyzqr_file >> x >> y >> z >> charge >> rad) {
		// ignore rest of the line...
		xyzqr_file.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

		// skip zero charges
		if (skip_zeroes && charge == 0.0) continue;

		net_charge += charge;

		Charge ch(Vector(x,y,z), charge, rad);

		//std::cout << " Charge magnitude " << charge << " at " << x << " " << y << " " << z << std::endl;
		charges.push_back(ch);
	}

	return net_charge;
}

#endif /* CHARGE_H_ */
