/*
 * base_vector.h
 *
 *  Created on: 21 Jul 2010
 *      Author: david
 */

#ifndef BASE_VECTOR_H_
#define BASE_VECTOR_H_

#define SQUARE_ROOT(x) sqrt(x);
#define R_SQUARE_ROOT(x) (1.0 / sqrt(x))
#define EXPONENTIAL exp

#include <ostream>
#include <sstream>
#include <cmath>

template <typename PrecisionType>
class __Vector
{
public:    
    
    PrecisionType x;
    PrecisionType y;
    PrecisionType z;

public:
    
    typedef __Vector<PrecisionType> thisType;
    
    
    __Vector() : x(0), y(0), z(0) {}
    __Vector(PrecisionType _x, PrecisionType _y, PrecisionType _z) : x(_x), y(_y), z(_z) {}

    __Vector(const __Vector<PrecisionType>& other) {
        x = other.x;
        y = other.y;
        z = other.z;
    }
    

    inline PrecisionType length() const
    {
        PrecisionType x2 = x * x;
        PrecisionType y2 = y * y;
        PrecisionType z2 = z * z;

        return SQUARE_ROOT(x2+y2+z2);

    };

    inline PrecisionType length2() const
    {
        PrecisionType x2 = x * x;
        PrecisionType y2 = y * y;
        PrecisionType z2 = z * z;

        return x2+y2+z2;

    };
    
    inline void normalise() { this->__divide(this->length()); }
    inline thisType normalised() const { 
        if (this->length2() == 0) { throw std::exception(); }
        return thisType(*this) / this->length(); 
    }

    template<typename other_precision>
    inline PrecisionType dot(const __Vector<other_precision>& v2) const
    {
        const __Vector<PrecisionType>& v1 = *this;
        return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
    }
    
    template<typename other_precision>
    inline thisType cross(const __Vector<other_precision>& v2) const 
    {
        const thisType& v1 = *this;
        
        __Vector<PrecisionType> result;
        result.x = v1.y * v2.z - v1.z * v2.y;
        result.y = v1.z * v2.x - v1.x * v2.z;
        result.z = v1.x * v2.y - v1.y * v2.x;

        return result;
    }

    inline PrecisionType operator()(int idx) const
    {
        switch (idx)
        {
        case 0:
            return x;
            break; // redundant
        case 1:
            return y;
            break; // redundant
        case 2:
            return z;
            break; // redundant
        default:
            throw std::exception();
            break; // redundant
        }
        return 0; // also redundant
    }

    inline PrecisionType& operator()(int idx) 
    {
        switch (idx)
        {
        case 0:
            return x;
            break; // redundant
        case 1:
            return y;
            break; // redundant
        case 2:
            return z;
            break; // redundant
        default:
            throw std::exception();
            break; // redundant
        }
        return x; // also redundant
    }

    inline thisType operator+(const thisType& other) const
    {
        return thisType(add(*this, other));
    }

    inline thisType operator-(const thisType& other) const
    {
        return thisType(subtract(*this, other));
    }

    inline thisType& operator*=(PrecisionType rhs)
    {
        x *= rhs;
        y *= rhs;
        z *= rhs;
        return *this;
    }

    inline thisType& operator/=(PrecisionType divisor)
    {
        x /= divisor;
        y /= divisor;
        z /= divisor;
        return *this;
    }

    template<typename other_precision>
    inline thisType& operator+=(const __Vector<other_precision>& rhs)
    {
        x += static_cast<PrecisionType>(rhs.x);
        y += static_cast<PrecisionType>(rhs.y);
        z += static_cast<PrecisionType>(rhs.z);
        return *this;
    }

    inline thisType& operator-=(const thisType &rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
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
        return thisType(x*multiplier, y*multiplier, z*multiplier);
    }

    template<typename other_precision>
    inline thisType& operator=(const __Vector<other_precision>& other)
    {
        x = static_cast<PrecisionType>(other.x);
        y = static_cast<PrecisionType>(other.y);
        z = static_cast<PrecisionType>(other.z);
        return *this;
    }

    inline thisType& operator=(PrecisionType val)
    {
        x = val;
        y = val;
        z = val;
        return *this;
    }
    
    inline void reset()
    {
        x = 0;
        y = 0;
        z = 0;
    }

    inline bool operator==(const thisType& other) const
    {
        return vectors_are_equal(*this, other);
    }

    inline bool operator!=(const thisType& other) const
    {
        return !vectors_are_equal(*this, other);
    }

    template<typename T>
    inline void toSpherical(T &rho, T &theta, T &phi) const
    {
        // converter from cartesian to spherical coordinates
        //PrecisionType x,y,z,rho,theta,phi
        //PrecisionType precision halfpi
        static const PrecisionType halfpi = 1.5707963267948966;

        rho = sqrt(x*x + y*y + z*z);
        theta = atan2(y,x);
        if (rho == 0.0){
            if (z < 0){
                phi = -halfpi;
            }
            else {
                phi = halfpi;
            }
        }
        else {
            assert(z/rho <= 1.0);
            phi = acos(z/rho);
        }
        return;
    }

    std::string __str__() const
    {
        std::ostringstream oss;
        oss << x << " " << y << " " << z << " ";
        return oss.str();
    }

protected:
    
    inline void __divide(PrecisionType divisor) {
        x /= divisor;
        y /= divisor;
        z /= divisor;
    }
    inline void __mult(PrecisionType multiple) {
        x *= multiple;
        y *= multiple;
        z *= multiple;
    }


};

template <typename PrecisionType, typename otherType>
inline __Vector<PrecisionType> operator*(otherType lhs, const __Vector<PrecisionType>& rhs)
{
    __Vector<PrecisionType> tmp(rhs);
    tmp *= lhs;
    return tmp;
}

template<typename PrecisionType>
inline __Vector<PrecisionType> operator-(const __Vector<PrecisionType>& other)
{
    return __Vector<PrecisionType>(-other.x, -other.y, -other.z);
}


template<typename PrecisionType> inline __Vector<PrecisionType> add(const __Vector<PrecisionType>& a, const __Vector<PrecisionType>& b)
{
    // Return a+b as __Vector<PrecisionType>
    __Vector<PrecisionType> a_plus_b;
    a_plus_b.x = a.x + b.x;
    a_plus_b.y = a.y + b.y;
    a_plus_b.z = a.z + b.z;
    return a_plus_b;
}

template<typename PrecisionType> inline __Vector<PrecisionType> subtract(const __Vector<PrecisionType>& a, const __Vector<PrecisionType>& b)
{
    // Return a-b as __Vector<PrecisionType>
    __Vector<PrecisionType> a_minus_b;
    a_minus_b.x = a.x - b.x;
    a_minus_b.y = a.y - b.y;
    a_minus_b.z = a.z - b.z;
    return a_minus_b;
}

template<typename PrecisionType> inline PrecisionType distance(const __Vector<PrecisionType>& p1, const __Vector<PrecisionType>& p2)
{
    PrecisionType dx = p1.x - p2.x;
    PrecisionType dy = p1.y - p2.y;
    PrecisionType dz = p1.z - p2.z;

    return SQUARE_ROOT(dx*dx + dy*dy + dz*dz);
}

template<typename PrecisionType> inline PrecisionType distance2(const __Vector<PrecisionType>& p1, const __Vector<PrecisionType>& p2)
{
    PrecisionType dx = p1.x - p2.x;
    PrecisionType dy = p1.y - p2.y;
    PrecisionType dz = p1.z - p2.z;

    return dx*dx + dy*dy + dz*dz;
}

template<typename PrecisionType> inline PrecisionType cosTheta(const __Vector<PrecisionType>& p1, const __Vector<PrecisionType>& p2)
{
    return p1.dot(p2) * R_SQUARE_ROOT(p1.length2()*p2.length2());
}

template<typename PrecisionType> inline bool vectors_are_equal(const __Vector<PrecisionType>& v1, const __Vector<PrecisionType>& v2)
{
	return (v1.x == v2.x && v1.y == v2.y && v1.z == v2.z);
}

// stream operators
template<typename PrecisionType>
inline std::ostream& operator<<(std::ostream& os, const __Vector<PrecisionType>& v)
{
    os << "(" << v.x << "," << v.y << "," << v.z << ")";
    return os;
}

#endif /* BASE_VECTOR_H_ */
