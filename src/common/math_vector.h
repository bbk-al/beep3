/*
 * math_vector.h
 *
 *  Created on: 21 Jul 2010
 *      Author: david
 */

#ifndef MATH_VECTOR_H_
#define MATH_VECTOR_H_

#ifdef PYTHON_MODULE
    #include <boost/python.hpp>
#endif

#include <iostream>
#include <ostream>
#include <sstream>
#include <cassert>
#include <limits>

#include "base_vector.h"

#ifdef __CHARMC__
#include <pup_stl.h>
#endif

const double pi = 3.141592653589793;
const double ONE_OVER_4PI = 1.0 / (4.0*pi);
const double PRECIS = std::numeric_limits<double>::epsilon();

class Quaternion
{
public:

    Quaternion() : a(1.0),b(0),c(0),d(0) {}
    Quaternion(double _a, double _b, double _c, double _d) :
        a(_a),
        b(_b),
        c(_c),
        d(_d) {}
    Quaternion(const Quaternion& other) : a(other.a), b(other.b), c(other.c), d(other.d) {}
    ~Quaternion() {}

    double a; // real component
    double b; // i coeff
    double c; // j coeff
    double d; // k coeff

    inline Quaternion conjugate() const { return Quaternion(a,-b,-c,-d); }
    inline Quaternion inverse() const { 
        return conjugate() / norm2(); 
    }
    
    inline Quaternion operator/(double divisor) const
    {
        return Quaternion(a/divisor, b/divisor, c/divisor, d/divisor);
    }

    template<typename T>
    inline void operator()(T &v1, T &v2, T &v3) const
    {
        T t2 =   a*b;
        T t3 =   a*c;
        T t4 =   a*d;
        T t5 =  -b*b;
        T t6 =   b*c;
        T t7 =   b*d;
        T t8 =  -c*c;
        T t9 =   c*d;
        T t10 = -d*d;

        T new_v1 = 2.0*( (t8 + t10)*v1 + (t6 -  t4)*v2 + (t3 + t7)*v3 ) + v1;
        T new_v2 = 2.0*( (t4 +  t6)*v1 + (t5 + t10)*v2 + (t9 - t2)*v3 ) + v2;
        T new_v3 = 2.0*( (t7 -  t3)*v1 + (t2 +  t9)*v2 + (t5 + t8)*v3 ) + v3;

        v1 = new_v1;
        v2 = new_v2;
        v3 = new_v3;

        return;
    }
    
    inline double norm2() const { return a*a + b*b + c*c + d*d; }
    inline double norm() const { return sqrt(a*a + b*b + c*c + d*d); }
    inline void normalise() { 
        double divisor = norm();
        a /= divisor;
        b /= divisor;
        c /= divisor;
        d /= divisor;
    }

#ifdef __CHARMC__
    void pup(PUP::er &p) {

        p | a;
        p | b;
        p | c;
        p | d;
    }
#endif

    inline std::string str() const
    {
    	std::ostringstream buf;
    	const Quaternion& q = *this;
    	buf << q.a << "," << q.b << "," << q.c << "," << q.d;
    	return buf.str();
    }

    //friend template<typename PrecisionType> std::ostream& operator<< (std::ostream &out, const Quaternion &q);

};

inline std::ostream& operator<<(std::ostream& os, const Quaternion& q)
{
    os << "(" << q.str() << ")";
    return os;
}

template <typename PrecisionType>
class VectorT : public __Vector<PrecisionType>
{

public:

    typedef VectorT<PrecisionType> thisType;
    typedef __Vector<PrecisionType> super;
    
    // possible downcast from double to float.
    // Wish we could template constructors... hmmm C++0x?
    template<typename another>
    operator VectorT<another>() const {
        VectorT<another> tmp(static_cast<PrecisionType>(super::x),
                             static_cast<PrecisionType>(super::y),
                             static_cast<PrecisionType>(super::z));
        return tmp;
    }

    VectorT() {}
    VectorT(const __Vector<float>& other) : super(other.x, other.y, other.z) {}
    VectorT(const __Vector<double>& other) : super(other.x, other.y, other.z) {}
    VectorT(PrecisionType _x, PrecisionType _y, PrecisionType _z){
        super::x = _x;
        super::y = _y;
        super::z = _z;
    }

    virtual ~VectorT() {}
    
public:

#ifdef __CHARMC__
    void pup(PUP::er &p) {

        p | super::x;
        p | super::y;
        p | super::z;
    }
#endif

    inline PrecisionType operator[](int idx) const
    {
        switch(idx)
        {
            case(0):
                return super::x;
            case(1):
                return super::y;
            case(2):
                return super::z;
            default:
                // should not get here
                throw;
        }
        // should not get here
        return super::x; // just to prevent compiler warnings...
    }

    inline PrecisionType dot(const thisType& v2) const
    {
        const __Vector<PrecisionType>& v1 = *this;
        return v1.dot(v2);
    }
    
    inline thisType cross(const thisType& v2) const
    {
        const __Vector<PrecisionType>& v1 = *this;
        return v1.cross(v2);
    }

    template<typename other_precision>
    inline thisType operator+(const __Vector<other_precision>& other) const
    {
        return thisType(super::x + other.x, 
                        super::y + other.y, 
                        super::z + other.z);
    }
    inline thisType __py_add(const thisType &other) const
    {
        return thisType(super::x + other.x, 
                        super::y + other.y, 
                        super::z + other.z);
    }

    template<typename other_precision>
    inline thisType operator-(const __Vector<other_precision>& other) const
    {
        return thisType(super::x - other.x, 
                        super::y - other.y, 
                        super::z - other.z);
    }
    inline thisType __py_subtract(const thisType &other) const
    {
        return thisType(super::x - other.x, 
                        super::y - other.y, 
                        super::z - other.z);
    }

    template<typename other_precision>
    inline thisType& operator+=(const __Vector<other_precision>& rhs)
    {
        super::x += static_cast<PrecisionType>(rhs.x);
        super::y += static_cast<PrecisionType>(rhs.y);
        super::z += static_cast<PrecisionType>(rhs.z);
        return *this;
    }

    inline thisType& operator-=(const thisType &rhs)
    {
        super::x -= rhs.x;
        super::y -= rhs.y;
        super::z -= rhs.z;
        return *this;
    }

    inline thisType operator/(PrecisionType divisor) const
    {
        thisType new_v(*this);
        new_v.__divide(divisor);
        return new_v;
    }

    inline thisType operator*(PrecisionType multiplier) const
    {
        return thisType(super::x*multiplier, super::y*multiplier, super::z*multiplier);
    }

    template<typename other_precision>
    inline thisType& operator=(const __Vector<other_precision>& other)
    {
        super::x = static_cast<PrecisionType>(other.x);
        super::y = static_cast<PrecisionType>(other.y);
        super::z = static_cast<PrecisionType>(other.z);
        return *this;    
    }

    inline thisType& operator=(double val)
    {
        super::x = static_cast<PrecisionType>(val);
        super::y = static_cast<PrecisionType>(val);
        super::z = static_cast<PrecisionType>(val);
        return *this;
    }

    inline std::string kinemage() const
    {
        std::ostringstream buf;
        buf << "{} P " << super::x << " " << super::y << " " << super::z << "\n";
        return buf.str();
    }

    // express this point as defined by the centre of rotation in the new frame; rot is a rotation operator from the old
    // coordinate frame into the new coordinate frame.
    void py_change_coordinate_frame(const thisType& centre_of_rotation, const Quaternion& rot, const thisType& centre_of_rotation_in_new_frame)
    {
        *this -= centre_of_rotation; // get relative to centre of rotation
        apply_rotation(rot);         // apply rotation to vector -- now points from centre of rotation in new frame
        *this += centre_of_rotation_in_new_frame;

    }

    // express this point as defined by the centre of rotation in the new frame; rot is a rotation operator from the old
    // coordinate frame into the new coordinate frame.
    void change_coordinate_frame(const thisType& centre_of_rotation, const Quaternion& rot, const thisType& centre_of_rotation_in_new_frame)
    {
        *this -= centre_of_rotation; // get relative to centre of rotation
        apply_rotation(rot);         // apply rotation to vector -- now points from centre of rotation in new frame
        *this += centre_of_rotation_in_new_frame;

    }

    // express this point as defined by the centre of rotation in the new frame; rot is a rotation operator from the old
    // coordinate frame into the new coordinate frame.
    thisType change_coordinate_frame(const thisType& centre_of_rotation, const Quaternion& rot, const thisType& centre_of_rotation_in_new_frame) const
    {
    	thisType new_v = *this;
        new_v -= centre_of_rotation; // get relative to centre of rotation
        new_v.apply_rotation(rot);         // apply rotation to vector -- now points from centre of rotation in new frame
        new_v += centre_of_rotation_in_new_frame;
        return new_v;

    }

    // for applying quaternions to vector in-situ
    void apply_rotation(const Quaternion& rot)
    {
        thisType& v = *this;
        rot(v.x, v.y, v.z);
        return;
    }

    // for applying quaternions to vectors
    thisType apply_rotation(const Quaternion& rot) const
    {
        thisType v(*this);
        v.apply_rotation(rot);
        return v;
    }

    // another alias for applyig a rotation- boost.python no good with overloaded funcs though
    inline void apply_quaternion(const Quaternion& rot) {
    	apply_rotation(rot);
    }

};

template<typename PrecisionType>
class PointNormalT
{
    
public:
    
    PointNormalT() : _pt(0,0,0), _n(0,0,0) {}
    PointNormalT(const VectorT<PrecisionType>& x, const VectorT<PrecisionType>& y) : _pt(x), _n(y) {}
    PointNormalT(const PointNormalT& other) : _pt(other._pt), _n(other._n) {}
    
    inline const VectorT<PrecisionType>& pt() const { return _pt; }
    inline const VectorT<PrecisionType>& n() const { return _n; }
    
    inline VectorT<PrecisionType>& pt() { return _pt; }
    inline VectorT<PrecisionType>& n() { return _n; }

    inline void change_coordinate_frame(const VectorT<PrecisionType>& centre, const Quaternion& rot, const VectorT<PrecisionType>& new_centre)
    {
        _pt.change_coordinate_frame(centre, rot, new_centre);
        _n.apply_rotation(rot);
    }

#ifdef __CHARMC__
    void pup(PUP::er &p) {

        p | _pt;
        p | _n;
    }
#endif

private:
    
    VectorT<PrecisionType> _pt;
    VectorT<PrecisionType> _n;
};

// default to double precision
typedef VectorT<double> Vector;
typedef PointNormalT<double> PointNormal;

class KahanVector : public Vector
{
    
public:

    KahanVector() : Vector(0,0,0), kahan(0,0,0) {}
    KahanVector(const KahanVector& other) : Vector(other), kahan(other.kahan) {}
    KahanVector(const Vector other) : Vector(other), kahan(0,0,0) {}

    inline KahanVector& operator+=(const KahanVector &rhs)
    {
        kahan += rhs.kahan;
        operator+=(static_cast<const Vector&>(rhs));
        return *this;
    }
    
    inline KahanVector& operator+=(const Vector &rhs)  //__attribute__ ((optimize(0)))
    {
        Vector& value = static_cast<Vector&>(*this);
        
        Vector y_vec(rhs - kahan);
        Vector t_vec(value + y_vec);
        kahan = (t_vec - value) - y_vec;
        value = t_vec;
        return *this;
    }
    
    inline KahanVector& operator-=(const Vector &rhs)
    {
        Vector other = rhs * -1;
        return this->operator+=(other);
    }
    
    inline KahanVector operator+(const KahanVector &rhs) const
    {
        KahanVector cpy(*this);
        cpy += rhs;
        return cpy;
    }

    inline KahanVector operator-(const KahanVector &rhs) const
    {
        KahanVector cpy(*this);
        cpy -= rhs;
        return cpy;
    }

    inline KahanVector operator*(double number) const
    {
        KahanVector cpy(*this);
        cpy *= number;
        return cpy;
    }

    inline KahanVector& operator*=(double rhs)
    {
        Vector::operator*=(rhs);
        kahan *= rhs;
        return *this;
    }
    
    inline KahanVector& operator/=(double &rhs)
    {
        Vector::operator/=(rhs);
        kahan /= rhs;
        return *this;
    }

    inline KahanVector& operator=(const KahanVector &other)
    {
        Vector::operator=(other);
        kahan = other.kahan;
        return *this;
    }

    inline KahanVector& operator=(const Vector &other)
    {
        Vector::operator=(other);
        kahan = Vector(0,0,0);
        return *this;
    }
    
    inline const Vector operator*() const { return *this; }
    inline Vector operator*() { return *this; }
    inline const Vector get_kahan() const { return kahan; }
    inline Vector get_kahan() { return kahan; }
    
protected:
    
    Vector kahan;
        
};

#endif /* MATH_VECTOR_H_ */
