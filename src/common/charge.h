/*
* charge.h
*
*  Created on: 21 Jul 2010
*      Author: david
*/

#ifndef CHARGE_H_
#define CHARGE_H_

#include "../common/math_vector.h"
#include <fstream>
#include <limits>
#include <string>
#include <vector>
#ifdef __CHARMC__
#include <pup_stl.h>
#endif

class Charge : public Vector
{
public:

    Charge() : Vector(0,0,0),  charge(0), radius(0) {}
    Charge(const Vector& v, double ch) : Vector(v), charge(ch), radius(0) {}
    Charge(const Vector& v, double ch, double r) : Vector(v), charge(ch), radius(r) {}
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
    ~Charge() {}

    inline Vector& position() { return static_cast<Vector&>(*this); }
    inline const Vector& position() const { return static_cast<const Vector&>(*this); }
    inline const Vector& py_get_position() const { return static_cast<const Vector&>(*this); }

    inline double get_charge() const { return charge; }
    inline double& get_charge() { return charge; }

    inline double get_charge(int) const { return charge; }
    inline double& get_charge(int) { return charge; }

    static double read_charges_from_file(const std::string& xyzqr_filename, std::vector<Charge>& charges, bool skip_zeroes=true)
    {
        // open the file
        std::ifstream xyzqr_file;
        xyzqr_file.open(xyzqr_filename.c_str());
        assert(xyzqr_file.good());

        double net_charge = 0.0;

        // read in xyzr data, and create charges
        double x,y,z,charge,rad;
        while(xyzqr_file >> x >> y >> z >> charge >> rad)
        {
            // ignore rest of the line...
            xyzqr_file.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

            // skip zero charges
            if (skip_zeroes && charge == 0.0) { continue; }

            net_charge += charge;

            Charge ch(Vector(x,y,z), charge, rad);

            //std::cout << " Charge magnitude " << charge << " at " << x << " " << y << " " << z << std::endl;
            charges.push_back(ch);
        }

        return net_charge;
    }

    double charge;
    double radius;

#ifdef __CHARMC__

    /// PUP Routine ///
    void pup(PUP::er &p) {

        Vector::pup(p);
        p | charge;
        p | radius;

    }
#endif

};

typedef Charge CharmChargeHolder;

#endif /* CHARGE_H_ */
