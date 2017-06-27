 /*      Author: david fallaize  Created on: 21 Jul 2010 */
 /*      Modified: adam light    on: 9 Mar 2012  */
/*! \file math_vector.h
 * \brief This module declares (and implements) vector and quaternion classes.
 *
 * Quaternions are used to implement rotations of Vectors, which in turn
 * are used to represent the positions of molecules (meshes and charges).
 * 
 *	This module contains the following public classes:
 *	- Quaternion -- representation of a + b i + c j + d k
 *	- VectorT -- templated representation of a 3-vector (x, y, z)
 *	- Vector -- defined as a VectorT of type double.
 *	- PointNormalT -- is a point vector and a normal vector (templated)
 *	- KahanVector -- is a vector with Kahan summation applied, which reduces
 *		errors in sums of a large number of terms.
 *	In general, these only implement operations actually used by BEEP and do
 *	not provide a general interface.
 *
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
#include <math.h>   // for M_PI
#include <random>

#include "base_vector.h"

#ifdef __CHARMC__
#include <pup_stl.h>
#endif

constexpr double pi = M_PI;
constexpr double ONE_OVER_4PI = 1.0 / (4.0*pi);
constexpr double PRECIS = std::numeric_limits<double>::epsilon();

//>>> There needs to be a math_vector.cpp to take away the silly inlines
// and this kludge...
#ifdef QUATERNION_REQUIRE_RANDOM
// These are needed to support randomisation but don't belong in a class
static std::random_device s_rd;		// obtain seed for random number engine
static std::mt19937_64 s_gen(s_rd());	// 64-bit mersenne_twister_engine
static std::uniform_real_distribution<double> s_dis(0.0, 1.0);	// Uniform [0,1)
#endif // QUATERNION_REQUIRE_RANDOM

class Quaternion;

template <typename PrecisionType>
class VectorT : public __Vector<PrecisionType> {

public:

    using thisType = VectorT<PrecisionType>;
    using super = __Vector<PrecisionType>;
    
	// No implicit special member functions in template class ...
    VectorT() = default;
    VectorT(const thisType&) = default;						//! Copy constructor
	constexpr thisType& operator=(const thisType&)
		= default;											//! Copy assignment
    VectorT(thisType&&) = default;							//! Move constructor
	constexpr thisType& operator=(thisType&&) = default;	//! Move assignment
    virtual ~VectorT() = default;

	//! Constructors from superclass
#ifdef DELETED
    VectorT(const __Vector<float>& other) : super(other.x, other.y, other.z) {}
    VectorT(const __Vector<double>& other) : super(other.x, other.y, other.z) {}
#else //  DELETED
	//! It is the responsibility of the caller to ensure that NewTypeT can
	//! sensibly cast to PrecisionType

	//! Copy from superclass
	template <typename NewTypeT>
    VectorT(const __Vector<NewTypeT>& other)
	: super(static_cast<PrecisionType>(other.x),
			static_cast<PrecisionType>(other.y),
			static_cast<PrecisionType>(other.z))
	{}
	//! Move from superclass
	template <typename NewTypeT>
    VectorT(__Vector<NewTypeT>&& other)
	: super(std::move(static_cast<PrecisionType>(other.x)),
			std::move(static_cast<PrecisionType>(other.y)),
			std::move(static_cast<PrecisionType>(other.z)))
	{}
#endif //  DELETED

	//! Assignment from other type
    template<typename NewTypeT>
    inline thisType& operator=(const __Vector<NewTypeT>& other);

    //! Direct Vector cast
#ifdef DELETED
    template<typename another>
    operator VectorT<another>() const {
        VectorT<another> tmp(static_cast<PrecisionType>(super::x),
                             static_cast<PrecisionType>(super::y),
                             static_cast<PrecisionType>(super::z));
        return tmp;
    }
#else // DELETED
	template <typename NewTypeT>
	inline operator VectorT<NewTypeT>();
#endif //  DELETED

	//! Constructor from coordinate values
    VectorT(PrecisionType _x, PrecisionType _y, PrecisionType _z) {
        super::x = _x;
        super::y = _y;
        super::z = _z;
    }
    
public:

#ifdef __CHARMC__
    void pup(PUP::er &p) {

        p | super::x;
        p | super::y;
        p | super::z;
    }
#endif

	//! Index operator to retrieve x as 0, y as 1, z as 2
    inline PrecisionType operator[](int idx) const;

	//! Dot and cross operations
    PrecisionType dot(const thisType& v2) const {
        const __Vector<PrecisionType>& v1 = *this;
        return v1.dot(v2);
    }
    
    thisType cross(const thisType& v2) const {
        const __Vector<PrecisionType>& v1 = *this;
        return v1.cross(v2);
    }

	//! Addition and subtraction of other types
	//! Responsibility of caller to ensure the casts are sensible.
    template<typename NewTypeT>
    thisType operator+(const __Vector<NewTypeT>& other) const {
        return thisType(super::x + other.x, 
                        super::y + other.y, 
                        super::z + other.z);
    }
#ifdef DELETED
    thisType __py_add(const thisType &other) const {
        return thisType(super::x + other.x, 
                        super::y + other.y, 
                        super::z + other.z);
    }
#endif //  DELETED

    template<typename NewTypeT>
    thisType operator-(const __Vector<NewTypeT>& other) const {
        return thisType(super::x - other.x, 
                        super::y - other.y, 
                        super::z - other.z);
    }
#ifdef DELETED
    thisType __py_subtract(const thisType &other) const {
        return thisType(super::x - other.x, 
                        super::y - other.y, 
                        super::z - other.z);
    }
#endif //  DELETED

	//! Addition and subtraction by other types
	//! Responsibility of caller to ensure the casts are sensible.
    template<typename NewTypeT>
    inline thisType& operator+=(const __Vector<NewTypeT>& rhs);

    template<typename NewTypeT>
    inline thisType& operator-=(const __Vector<NewTypeT>& rhs);

	//! Multiplication and division by scalars
    inline thisType operator/(PrecisionType divisor) const;

    thisType operator*(PrecisionType multiplier) const {
        return thisType(super::x*multiplier,
						super::y*multiplier,
						super::z*multiplier);
    }

    inline thisType& operator=(double val);

    inline std::string kinemage() const;

#ifdef DELETED  // Now defined in pybeep.cpp
    void py_change_coordinate_frame(const thisType& centre_of_rotation, const Quaternion& rot, const thisType& centre_of_rotation_in_new_frame)
	{
		change_coordinate_frame(centre_of_rotation, rot, centre_of_rotation_in_new_frame);
	}
#endif // DELETED

    // express this point as defined by the centre of rotation in the new 
    // frame; rot is a rotation operator from the old coordinate frame into
    // the new coordinate frame.
    inline void change_coordinate_frame(
		const thisType& centre_of_rotation,
		const Quaternion& rot,
		const thisType& centre_of_rotation_in_new_frame);
    // const version:
    inline thisType change_coordinate_frame(
		const thisType& centre_of_rotation,
		const Quaternion& rot,
		const thisType& centre_of_rotation_in_new_frame) const;

    // for applying quaternions to vector in-situ
    inline void apply_rotation(const Quaternion& rot);
    // for applying quaternions to vectors
    inline thisType apply_rotation(const Quaternion& rot) const;

#ifdef DELETED  // Now defined in pybeep.cpp
    // another alias for applyig a rotation- boost.python no good with overloaded funcs though
    inline void apply_quaternion(const Quaternion& rot) { apply_rotation(rot); }
#endif  //  DELETED

};

template<typename PrecisionType>
class PointNormalT
{
public:
    
    PointNormalT() = default;
    PointNormalT(const PointNormalT<PrecisionType>&)
		= default;											//! Copy constructor
	constexpr PointNormalT<PrecisionType>& operator=
		(const PointNormalT<PrecisionType>&) = default;		//! Copy assignment
    PointNormalT(PointNormalT<PrecisionType>&&) = default;	//! Move constructor
	constexpr PointNormalT<PrecisionType>& operator=
		(PointNormalT<PrecisionType>&&) = default;			//! Move assignment

	//! Constructor from point and normal vectors
    PointNormalT(const VectorT<PrecisionType>& x,
				 const VectorT<PrecisionType>& y) : _pt(x), _n(y) {}
    
	// accessors
    const VectorT<PrecisionType>& pt() const { return _pt; }
    const VectorT<PrecisionType>& n() const { return _n; }
    
	// modify accessors
    VectorT<PrecisionType>& pt() { return _pt; }
    VectorT<PrecisionType>& n() { return _n; }

    inline void change_coordinate_frame(
		const VectorT<PrecisionType>& centre,
		const Quaternion& rot,
		const VectorT<PrecisionType>& new_centre);

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


class Quaternion {
public:

	// Directly accessible coefficients:
    double a; // real component
    double b; // i coeff
    double c; // j coeff
    double d; // k coeff

	//! Default constructor to *unit* quaternion
    Quaternion() : a(1.0),b(0),c(0),d(0) {}

	//! Constructor from coefficients
    Quaternion(double _a, double _b, double _c, double _d)
		: a(_a), b(_b), c(_c), d(_d)
	{}

    Quaternion(const Quaternion&) = default;			//! Copy constructor
	Quaternion& operator=(const Quaternion&) = default;	//! Copy assignment
    Quaternion(Quaternion&& other) = default;			//! Move constructor
	Quaternion& operator=(Quaternion&&) = default;		//! Move assignment

	//! Destructor
    virtual ~Quaternion() = default;

	// Useful operations
    bool operator==(const Quaternion& q) const { return equals(*this, q); }
    bool operator!=(const Quaternion& q) const { return !equals(*this, q); }
	bool equals(const Quaternion& q1, const Quaternion& q2) const {
		return (q1.a == q2.a && q1.b == q2.b && q1.c == q2.c && q1.d == q2.d);
	}
    Quaternion conjugate() const { return Quaternion(a,-b,-c,-d); }
    Quaternion inverse() const { return conjugate() / norm2(); }
    
    Quaternion operator/(double divisor) const {
        return Quaternion(a/divisor, b/divisor, c/divisor, d/divisor);
    }

	// Quaternion product
    inline Quaternion operator*(const Quaternion& q) const;
    inline Quaternion& operator*=(const Quaternion& q);

	// Function operators:  Quaternion q(...);  q(v1, v2, v3);  q(v);
	// These apply the quaternion rotation to the vector provided
    template<typename T>
    inline void operator()(T &v1, T &v2, T &v3) const;

    template<typename T>
    void operator()(VectorT<T> &v) const
    {
		return T(v.x, v.y, v.z);
	}
    
    double norm2() const { return a*a + b*b + c*c + d*d; }
    double norm() const { return sqrt(a*a + b*b + c*c + d*d); }
    inline void normalise();

#ifdef __CHARMC__
    void pup(PUP::er &p) {

        p | a;
        p | b;
        p | c;
        p | d;
    }
#endif

    inline std::string str() const;

#ifdef QUATERNION_REQUIRE_RANDOM
	// Generate a random rotation with a uniform distribution
	static inline Quaternion rand();
#endif // QUATERNION_REQUIRE_RANDOM

    //friend template<typename PrecisionType> std::ostream& operator<< (std::ostream &out, const Quaternion &q);

};

inline std::ostream& operator<<(std::ostream& os, const Quaternion& q) {
    os << "(" << q.str() << ")";
    return os;
}


class KahanVector : public Vector {
    
public:

    KahanVector() = default;
    KahanVector(const KahanVector&) = default;				//! Copy constructor
	KahanVector& operator=(const KahanVector&) = default;	//! Copy assignment
    KahanVector(KahanVector&& other) = default;				//! Move constructor
	KahanVector& operator=(KahanVector&&) = default;		//! Move assignment
	~KahanVector() = default;

	//! Constructor from vector
    KahanVector(const Vector& other) : Vector(other), kahan(0,0,0) {}
    KahanVector(Vector&& other) : Vector(other), kahan(0,0,0) {}

	//! Assignment from vector
    inline KahanVector& operator=(const Vector &other);

	// Operators
    inline KahanVector operator+(const KahanVector &rhs) const;
    inline KahanVector& operator+=(const KahanVector &rhs);
    inline KahanVector& operator+=(const Vector &rhs);  //__attribute__ ((optimize(0)))
    inline KahanVector operator-(const KahanVector &rhs) const;
    inline KahanVector& operator-=(const Vector &rhs);
    inline KahanVector operator*(double rhs) const;
    inline KahanVector& operator*=(double rhs);
    inline KahanVector operator/(double rhs) const;
    inline KahanVector& operator/=(double rhs);

	// ContentsOf as Vector
    const Vector operator*() const { return *this; }
    Vector operator*() { return *this; }

	// get_ methods
    const Vector& get_kahan() const { return kahan; }
    Vector& get_kahan() { return kahan; }
    
protected:
    
    Vector kahan;
        
};


// Inlined methods

// VectorT
template<class PrecisionType>
template<class NewTypeT>
inline VectorT<PrecisionType>::operator VectorT<NewTypeT>()
{
	VectorT<NewTypeT> tmp;
	tmp.x = static_cast<PrecisionType>(super::x);
	tmp.y = static_cast<PrecisionType>(super::y);
	tmp.z = static_cast<PrecisionType>(super::z);
	return tmp;
}

template<class PrecisionType>
inline PrecisionType VectorT<PrecisionType>::operator[](int idx) const
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

template<class PrecisionType>
template<typename NewTypeT>
inline VectorT<PrecisionType>& VectorT<PrecisionType>::operator+=
	(const __Vector<NewTypeT>& rhs)
{
	super::x += static_cast<PrecisionType>(rhs.x);
	super::y += static_cast<PrecisionType>(rhs.y);
	super::z += static_cast<PrecisionType>(rhs.z);
	return *this;
}

template<class PrecisionType>
template<typename NewTypeT>
inline VectorT<PrecisionType>& VectorT<PrecisionType>::operator-=
	(const __Vector<NewTypeT>& rhs)
{
	super::x -= static_cast<PrecisionType>(rhs.x);
	super::y -= static_cast<PrecisionType>(rhs.y);
	super::z -= static_cast<PrecisionType>(rhs.z);
	return *this;
}

template<class PrecisionType>
inline VectorT<PrecisionType> VectorT<PrecisionType>::operator/
	(PrecisionType divisor) const
{
	thisType new_v(*this);
	new_v.__divide(divisor);
	return new_v;
}

template<class PrecisionType>
template<class NewTypeT>
inline VectorT<PrecisionType>& VectorT<PrecisionType>::operator=
	(const __Vector<NewTypeT>& other)
{
	super::x = static_cast<PrecisionType>(other.x);
	super::y = static_cast<PrecisionType>(other.y);
	super::z = static_cast<PrecisionType>(other.z);
	return *this;    
}

template<class PrecisionType>
inline VectorT<PrecisionType>& VectorT<PrecisionType>::operator=(double val)
{
	super::x = static_cast<PrecisionType>(val);
	super::y = static_cast<PrecisionType>(val);
	super::z = static_cast<PrecisionType>(val);
	return *this;
}

template<class PrecisionType>
inline std::string VectorT<PrecisionType>::kinemage() const
{
	std::ostringstream buf;
	buf << "{} P "
		<< super::x << " " << super::y << " " << super::z << "\n";
	return buf.str();
}

template<class PrecisionType>
inline void VectorT<PrecisionType>::change_coordinate_frame(
	const thisType& centre_of_rotation,
	const Quaternion& rot,
	const thisType& centre_of_rotation_in_new_frame)
{
	*this -= centre_of_rotation; // get relative to centre of rotation
	apply_rotation(rot);         // apply rotation to vector
					// -- now points from centre of rotation in new frame
	*this += centre_of_rotation_in_new_frame;
}

template<class PrecisionType>
inline VectorT<PrecisionType> VectorT<PrecisionType>::change_coordinate_frame(
	const VectorT<PrecisionType>& centre_of_rotation,
	const Quaternion& rot,
	const VectorT<PrecisionType>& centre_of_rotation_in_new_frame) const
{
	thisType new_v = *this;
#ifndef __DELETED_
	new_v.change_coordinate_fram(centre_of_rotation, rot,
								centre_of_rotation_in_new_frame);
#else //  __DELETED_
	new_v -= centre_of_rotation; // get relative to centre of rotation
	new_v.apply_rotation(rot);         // apply rotation to vector -- now points from centre of rotation in new frame
	new_v += centre_of_rotation_in_new_frame;
#endif //  __DELETED_
	return new_v;
}

template<class PrecisionType>
inline void VectorT<PrecisionType>::apply_rotation(const Quaternion& rot)
{
	thisType& v = *this;
	rot(v.x, v.y, v.z);
	return;
}

template<class PrecisionType>
inline VectorT<PrecisionType> VectorT<PrecisionType>::apply_rotation
	(const Quaternion& rot) const
{
	thisType v(*this);
	v.apply_rotation(rot);
	return v;
}

// PointNormalT
template<class PrecisionType>
inline void PointNormalT<PrecisionType>::change_coordinate_frame(
	const VectorT<PrecisionType>& centre,
	const Quaternion& rot,
	const VectorT<PrecisionType>& new_centre)
{
	_pt.change_coordinate_frame(centre, rot, new_centre);
	_n.apply_rotation(rot);
}

// Quaternions
template<typename T>
inline void Quaternion::operator()(T &v1, T &v2, T &v3) const
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
}

inline void Quaternion::normalise() { 
	double divisor = norm();
	a /= divisor;
	b /= divisor;
	c /= divisor;
	d /= divisor;
}

inline std::string Quaternion::str() const
{
	std::ostringstream buf;
	const Quaternion& q = *this;
	buf << q.a << "," << q.b << "," << q.c << "," << q.d;
	return buf.str();
}

inline Quaternion Quaternion::operator*(const Quaternion& q) const {
	Quaternion p(*this);
	p.a = p.a*q.a - p.b*q.b - p.c*q.c - p.d*q.d;
	p.b = p.a*q.b + p.b*q.a + p.c*q.d - p.d*q.c;
	p.c = p.a*q.c + p.c*q.a + p.d*q.b - p.b*q.d;
	p.d = p.a*q.d + p.d*q.a + p.b*q.c - p.c*q.b;
	return p;
}

inline Quaternion& Quaternion::operator*=(const Quaternion& q) {
	*this = *this * q;
	return *this;
}
 
#ifdef QUATERNION_REQUIRE_RANDOM
inline Quaternion Quaternion::rand() {
	// Using global distribution and generator specified after includes
	double s = s_dis(s_gen);  // squared sine at this point
	double c = sqrt(1 - s);
	s = sqrt(s);
	double a = s_dis(s_gen) * 2*M_PI;
	double b = s_dis(s_gen) * 2*M_PI;
	return Quaternion(cos(a) * s, sin(b) * c, cos(b) * c, sin(a) * s);
}
#endif // QUATERNION_REQUIRE_RANDOM

// KahanVector
inline KahanVector& KahanVector::operator=(const Vector &other) {
	Vector::operator=(other);
	kahan = Vector(0,0,0);
	return *this;
}

inline KahanVector KahanVector::operator+(const KahanVector &rhs) const
{
	KahanVector cpy(*this);
	cpy += rhs;
	return cpy;
}

inline KahanVector& KahanVector::operator+=(const KahanVector &rhs)
{
	kahan += rhs.kahan;
	operator+=(static_cast<const Vector&>(rhs));
	return *this;
}

inline KahanVector& KahanVector::operator+=(const Vector &rhs)  //__attribute__ ((optimize(0)))
{
	Vector& value = static_cast<Vector&>(*this);
	
	Vector y_vec(rhs - kahan);
	Vector t_vec(value + y_vec);
	kahan = (t_vec - value) - y_vec;
	value = t_vec;
	return *this;
}

inline KahanVector KahanVector::operator-(const KahanVector &rhs) const
{
	KahanVector cpy(*this);
	cpy -= rhs;
	return cpy;
}

inline KahanVector& KahanVector::operator-=(const Vector &rhs)
{
	Vector other = rhs * -1;
	return this->operator+=(other);
}

inline KahanVector KahanVector::operator*(double rhs) const {
	KahanVector cpy(*this);
	cpy *= rhs;
	return cpy;
}

inline KahanVector& KahanVector::operator*=(double rhs) {
	Vector::operator*=(rhs);
	kahan *= rhs;
	return *this;
}

inline KahanVector KahanVector::operator/(double rhs) const {
	KahanVector cpy(*this);
	cpy /= rhs;
	return cpy;
}

inline KahanVector& KahanVector::operator/=(double rhs) {
	Vector::operator/=(rhs);
	kahan /= rhs;
	return *this;
}

#endif /* MATH_VECTOR_H_ */
