/*
* quad_point.h
*
*  Created on: 21 Jul 2010
*      Author: david
*/

#ifndef QUAD_POINT_H_
#define QUAD_POINT_H_

#include "../common/math_vector.h"
#include <boost/math/special_functions/fpclassify.hpp>
#include <vector>
#include <sstream>
#include <string>
#include <cstring>
#include <ostream>
#include <exception>

#ifdef __CHARMC__
#include <charm++.h>
#include <pup_stl.h>
#endif


// holder for the abcissa/weights of actual quadrature rules
class QuadratureRule
{
public:

    QuadratureRule(unsigned int N, double x[], double y[], double w[])
    {
        num_points = N;
        quad_points_x = new double[N];
        quad_points_y = new double[N];
        quad_weights = new double[N];

        double total_weight = 0;
        for (unsigned int ii=0; ii < N; ++ii)
        {
            quad_points_x[ii] = x[ii];
            quad_points_y[ii] = y[ii];
            quad_weights[ii] = w[ii];
            total_weight += w[ii];
        }
        
        //std::cout << "Total quadrature weight: " << total_weight << std::endl;
        for (unsigned int ii=0; ii < N; ++ii)
        {
            quad_weights[ii] /= total_weight;
        }
        
    }

    ~QuadratureRule() {
        delete[] quad_points_x;
        delete[] quad_points_y;
        delete[] quad_weights;
    }

    unsigned int num_points;
    double *quad_points_x;
    double *quad_points_y;
    double *quad_weights;
};

class BadQuadPoint : public std::exception 
{
public:
    BadQuadPoint() : std::exception() {}
    inline std::string str() const {
        return "Someone tried to create a QuadPoint with NaN's.";
    }
};

inline std::ostream& operator<<(std::ostream& os, const BadQuadPoint& bad_qp_exception)
{
    os << bad_qp_exception.str() << std::endl;
    return os;
}

template<typename PrecisionType>
class QuadPointT
{
    // number of pieces of data stored in quad point
    static const unsigned int SIZE = 8;
    
public:

    static const unsigned int num_bytes = sizeof(PrecisionType)*SIZE;
    
    typedef QuadPointT<PrecisionType> thisType;
    typedef VectorT<PrecisionType> VecType;
    typedef __Vector<PrecisionType> VecStructType;

    QuadPointT() {
        assert(sizeof(VecStructType) == 3*sizeof(PrecisionType));
    }
    ~QuadPointT() {}
    QuadPointT(const VecType& _pt, const VecType& _normal, PrecisionType _weight)
    {
        for (int xx=0; xx < SIZE; ++xx)
        {
            flat_data[xx]=0;
        }
        assert(sizeof(VecStructType) == 3*sizeof(PrecisionType));
        pt() = _pt;
        normal() = _normal;
        weight() = _weight;
        
        for (int ii=0; ii < SIZE; ++ii)
        {
            if (boost::math::isnan(flat_data[ii]))
            {
                //std::cerr << "Error: Tried to construct a QuadPoint with NaNs..." << std::endl;
                throw BadQuadPoint();
            }
        }
        
    }
    
    // this is deprecated
    QuadPointT(const std::string& ascii_serialized) {
        
        assert(sizeof(VecStructType) == 3*sizeof(PrecisionType));
        
        // create quadrature point from line in serialized file
        std::stringstream stream_buffer(ascii_serialized);
        VecType& p = pt();
        VecType& n = normal();
        PrecisionType& w = weight();
        
        stream_buffer >> p.x >> p.y >> p.z >> n.x >> n.y >> n.z >> w;
        
    }

    // copy constructor
    QuadPointT(const thisType& other)
    {
        std::memcpy(flat_data, other.flat_data, sizeof(PrecisionType)*SIZE);
    }

    // copy construct with rotation to a new coordinate frame
    QuadPointT(const thisType& other, const VecType& centre_of_rotation, const Quaternion& rotation, const VecType& xyz_offset)
    {
        *this = other;
        assert(sizeof(VecStructType) == 3*sizeof(PrecisionType));
        change_coordinate_frame(centre_of_rotation, rotation, xyz_offset);
    }

    // polymorphism (sortof) with Vector types
    operator VectorT<PrecisionType>& () { return pt(); }
    operator const VectorT<PrecisionType>& () const { return pt(); }

    inline void change_coordinate_frame(const Vector& centre_of_rotation_old_frame,
                                        const Quaternion& rot_from_old_to_new,
                                        const Vector& centre_of_rotation_new_frame)
    {
        Vector _pt(pt());
        Vector _n(normal());
        
        _pt.change_coordinate_frame(centre_of_rotation_old_frame, rot_from_old_to_new, centre_of_rotation_new_frame);
        _n.apply_rotation(rot_from_old_to_new); // normal vector just gets rotated
        
        pt() = _pt;
        normal() = _n;
        
        return;
    }

    inline thisType& operator=(const VecType& other) {
        VecType& _pt = pt();
        _pt = other;
        return *this;
    }

    inline thisType& operator=(const thisType& other) {
        memcpy(flat_data, other.flat_data, sizeof(PrecisionType)*SIZE);
        return *this;
    }

    // evil accessor functions!  (These are like this so that I can directly
    // memcpy from QuadPoint objects to OpenCL quadpoint structs. So I need
    // a nice chunk of contiguous flat memory.)

    inline VecStructType& pt() { return *reinterpret_cast<VecStructType*>(&(flat_data[0])); }
    inline VecStructType& normal() { return *reinterpret_cast<VecStructType*>(&(flat_data[3])); }
    inline PrecisionType& weight() { return flat_data[6]; }

    inline const VecStructType& pt() const { return *reinterpret_cast<const VecStructType*>(&(flat_data[0])); }
    inline const VecStructType& normal() const { return *reinterpret_cast<const VecStructType*>(&(flat_data[3])); }
    inline const PrecisionType& weight() const { return flat_data[6]; }

    // python accessors (doesn't like overloading, tsk).
    //inline const VecStructType& __py_pt() const { return *reinterpret_cast<const VecStructType*>(&(flat_data[0])); }
    //inline const VecStructType& __py_normal() const { return *reinterpret_cast<const VecStructType*>(&(flat_data[3])); }
    //inline const PrecisionType& __py_weight() const { return flat_data[6]; }

    std::string str() const
    {
        std::stringstream buf;
        buf << pt() << " "
            << normal() << " "
            << weight();
        return buf.str();

    }

#ifdef __CHARMC__
    void pup(PUP::er &p) {
        p(flat_data,SIZE);
    }
#endif

    const void* bytes() const { return flat_data; }

private:
    
    PrecisionType flat_data[SIZE];

};

template<typename PrecisionType>
inline std::ostream& operator<<(std::ostream& os, const QuadPointT<PrecisionType>& qp)
{
    os << "(" << qp.str() << ")";
    return os;
}

// default to single-precision quadrature points
typedef QuadPointT<float> QuadPoint;

// typedef a list of QuadPoints (float precision)
typedef std::vector<QuadPoint> QuadList;

#endif /* QUAD_POINT_H_ */
