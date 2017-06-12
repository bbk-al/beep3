/*
* multipole_holder.h
*
*  Created on: 16 Aug 2010
*      Author: david
*/

#ifndef MULTIPOLE_HOLDER_H_
#define MULTIPOLE_HOLDER_H_

#include <complex>
#include <boost/shared_array.hpp>

// TODO: replace these with Charm asserts (which destroy objects properly)
#include <cassert>
#include <sstream>
#include <iostream>
#include <string>
#ifdef __CHARMC__
#include <pup_stl.h>
#endif

static const unsigned short STRIDE_ARRAY[32] = {0, 1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66, 78, 91, 105, 120, 136, 153, 171, 190, 210, 231, 253, 276, 300, 325, 351, 378, 406, 435, 465, 496};

template<unsigned short size, typename vartype=double>
class TriangularMemory
{
public:

    static const unsigned short SIZE = size;
    static const unsigned short HOLDER_SIZE = ((size+1)*(size+2))/2;
    TriangularMemory() { reset(); }

    inline vartype& operator()(short n, short m)
    {
#ifdef LTL_RANGE_CHECKING
        assert(n <= size);
        assert(n >=0);
        assert(abs(m) <= n);
        assert(m >=0);
#endif
        return data[m+ STRIDE_ARRAY[n]];
    }

    inline const vartype& operator() (short n, short m) const
    {
#ifdef LTL_RANGE_CHECKING
        assert(n <= size);
        assert(n >=0);
        assert(abs(m) <= n);
        assert(m >=0);
#endif
        return data[m+ STRIDE_ARRAY[n]];
    }

    // re-inits all values to zero
    inline void reset()
    {
        memset(data,0,sizeof(data));
#ifdef USE_KAHAN
        memset(kahan,0,sizeof(kahan));
#endif
        return;
    }

    std::string str() const
    {
        std::ostringstream oss;
        for (unsigned short n=0; n <= size; ++n)
        {
            for (unsigned short m=0; m <= n; ++m)
            {
                const vartype& num = (*this)(n,m);
                oss << "(n=" << n << ",m=" << m << ") " << num << " ";
            }
            oss << "\n";
        }
        return oss.str();
    }

#ifdef __CHARMC__
    /// PUP Routine ///
    void pup(PUP::er &p) {

        // Array of std::complex<double>'s
        for (unsigned short i=0; i < HOLDER_SIZE; ++i)
        {
            p | data[i];
#ifdef USE_KAHAN
            p | kahan[i];
#endif
        }
    }
#endif

protected:

    // triangular memory to hold multipole expansion terms
    vartype data[HOLDER_SIZE];
#ifdef USE_KAHAN
    vartype kahan[HOLDER_SIZE];
#endif
};

// stream operator for TriangularMemory things
template<unsigned short size, typename vartype>
inline std::ostream& operator<<(std::ostream& os, const TriangularMemory<size,vartype>& triangular_memory)
{
    os << triangular_memory.str();
    return os;
}

template<int NUM_TERMS>
class BaseMultipoleHolder : public TriangularMemory< NUM_TERMS, std::complex<double> >
{

public:

    typedef TriangularMemory< NUM_TERMS, std::complex<double> > super;

    BaseMultipoleHolder() : super() {}

    BaseMultipoleHolder(const BaseMultipoleHolder<NUM_TERMS>& other) : super(other)
    {
        assert(sizeof(super::data) == sizeof(super::vartype)*super::HOLDER_SIZE);
        memcpy(super::data, other.data, sizeof(super::data));
#ifdef USE_KAHAN
        memcpy(kahan, other.kahan, sizeof(kahan));
#endif
    }
    inline BaseMultipoleHolder<NUM_TERMS>& operator= (const BaseMultipoleHolder<NUM_TERMS>& other)
    {
        if (this != &other)
        {
            memcpy(super::data, other.data, sizeof(super::data));
#ifdef USE_KAHAN
            memcpy(kahan, other.kahan, sizeof(kahan));
#endif
        }

        return *this;
    }

    inline BaseMultipoleHolder<NUM_TERMS>& operator= (std::complex<double> val)
    {
        for (size_t i=0, lim=sizeof(super::data)/sizeof(std::complex<double>); i < lim; ++i)
        {
            super::data[i] = val;
        }
#ifdef USE_KAHAN
        memset(kahan, 0, sizeof(kahan));
#endif
        return *this;
    }

    inline BaseMultipoleHolder<NUM_TERMS> operator*(double mult) const
    {
        BaseMultipoleHolder<NUM_TERMS> cpy(*this);
        for (size_t i=0, lim=sizeof(super::data)/sizeof(std::complex<double>); i < lim; ++i)
        {
            cpy.data[i] *= mult;
        }
        return cpy;
    }

    inline BaseMultipoleHolder<NUM_TERMS>& operator+=(const BaseMultipoleHolder<NUM_TERMS> &rhs)
    {
        for (unsigned short i=0; i < super::HOLDER_SIZE; ++i)
        {
#ifndef USE_KAHAN
            super::data[i] += rhs.data[i];
#else
            std::complex<double>& sum = data[i];
            std::complex<double> val = rhs.data[i];
            std::complex<double>& kh = kahan[i];
            kh += rhs.kahan[i];
            {
                // kahan summation (aka compensated addition)
                std::complex<double> tmp_y = val - kh;
                std::complex<double> tmp_t = sum + tmp_y;
                kh = (tmp_t - sum) - tmp_y;
                sum = tmp_t;
            }
#endif
        }

        return *this;
    }

    inline BaseMultipoleHolder<NUM_TERMS>& operator/=(double val)
    {
        for (unsigned short i=0; i < super::HOLDER_SIZE; ++i)
        {
            super::data[i] /= val;
#ifdef USE_KAHAN
            kahan[i] /= val;
#endif
        }
        return *this;
    }

    inline bool operator==(double val) const
    {
        const std::complex<double> cval(val,0);
        for (unsigned short i=0; i < super::HOLDER_SIZE; ++i)
        {
            if (super::data[i] != cval) { return false; }
        }
        return true;
    }

    inline bool operator!=(double val) const
    {
        return !(*this == val);
    }

};

template<int NUM_LAMBDAS, int NUM_WAVES>
class BasePlaneWaveHolder
{

public:

    static const unsigned short SIZE_PLANE_WAVE_VECTOR = NUM_WAVES;

    std::complex<double> data[SIZE_PLANE_WAVE_VECTOR];

    BasePlaneWaveHolder() {reset();}
    
    BasePlaneWaveHolder(const BasePlaneWaveHolder<NUM_LAMBDAS, NUM_WAVES>& other)
    {
        memcpy(data, other.data, sizeof(data));
    }

    inline BasePlaneWaveHolder<NUM_LAMBDAS, NUM_WAVES>& operator= (const BasePlaneWaveHolder<NUM_LAMBDAS, NUM_WAVES>& other)
    {
        if (this != &other)
        {
            memcpy(data, other.data, sizeof(data));
        }
        return *this;
    }

    BasePlaneWaveHolder<NUM_LAMBDAS, NUM_WAVES>& operator+= (const BasePlaneWaveHolder<NUM_LAMBDAS, NUM_WAVES>& other)
    {
        for (unsigned short i=0; i < SIZE_PLANE_WAVE_VECTOR; ++i)
        {
            this->data[i] += other.data[i];
        }
        return *this;
    }

    inline bool is_empty() const
    {
        const std::complex<double> zero(0,0);
        for (unsigned short i=0; i < SIZE_PLANE_WAVE_VECTOR; ++i)
        {
            if (data[i] != zero)
            {
                return false;
            }
        }
        return true;
    }

    inline void reset()
    {
        memset(data,0,sizeof(data));

    }

    std::complex<double>& operator() (unsigned int idx)
    {
#ifdef LTL_RANGE_CHECKING
        assert(idx < SIZE_PLANE_WAVE_VECTOR);
#endif
        return data[idx];
    }

    BasePlaneWaveHolder<NUM_LAMBDAS, NUM_WAVES>& operator=(const double &val)
    {
        *this = std::complex<double>(val,0);
        return *this;
    }

    BasePlaneWaveHolder<NUM_LAMBDAS, NUM_WAVES>& operator=(const std::complex<double> &val)
    {
        for (unsigned short ii=0; ii < SIZE_PLANE_WAVE_VECTOR; ++ii)
        {
            data[ii] = val;
        }
        return *this;
    }

    const std::complex<double>& operator() (unsigned int idx) const
    {
#ifdef LTL_RANGE_CHECKING
        assert(idx < SIZE_PLANE_WAVE_VECTOR);
#endif
        return data[idx];
    }

    std::string str() const
    {
        std::ostringstream oss;
        for (unsigned short ii=0; ii < SIZE_PLANE_WAVE_VECTOR; ++ii)
        {
            std::complex<double> num = data[ii];
            oss << "(i=" << ii << ") " << num << "\n";
        }
        return oss.str();
    }

#ifdef __CHARMC__
    /// PUP Routine ///
    void pup(PUP::er &p) {

        // Array of std::complex<double>'s
        for (size_t i=0; i < SIZE_PLANE_WAVE_VECTOR; ++i)
        {
            p | data[i];
        }
    }
#endif

};

template <int multiplicity, typename HType>
class MultiHolder
{
    HType holders[multiplicity];

public:

    inline void operator= (const MultiHolder<multiplicity,HType>& other)
    {
        for (unsigned short ii=0; ii < multiplicity; ++ii)
        {
            (*this)[ii] = other[ii];
        }
    }

    inline HType& operator[] (unsigned int idx)
    {
        assert(idx < multiplicity);
        return holders[idx];
    }

    inline const HType& operator[] (unsigned int idx) const
    {
        assert(idx < multiplicity);
        return holders[idx];
    }

    inline HType& operator() (unsigned int idx)
    {
        assert(idx < multiplicity);
        return holders[idx];
    }

    inline const HType& operator() (unsigned int idx) const
    {
        assert(idx < multiplicity);
        return holders[idx];
    }

    // re-inits all values to zero
    inline void reset()
    {
        for (unsigned short ii=0; ii < multiplicity; ++ii)
        {
            holders[ii].reset();
        }
        return;
    }

    inline MultiHolder<multiplicity,HType>& operator+=(const MultiHolder<multiplicity,HType> &rhs)
    {
        for (unsigned short ii=0; ii < multiplicity; ++ii)
        {
            holders[ii] += rhs[ii];
        }
        return *this;
    }

#ifdef __CHARMC__
    /// PUP Routine ///
    void pup(PUP::er &p) {

        // Array of std::complex<double>'s
        for (unsigned short ii=0; ii < multiplicity; ++ii)
        {
            p | holders[ii];
        }
    }
#endif

    std::string str() const
    {
        std::ostringstream oss;
        for (unsigned short ii=0; ii < multiplicity; ++ii)
        {
            oss << holders[ii].str();
        }
        return oss.str();
    }

};

// Template specialisation for single holders -- revert to base type
template <typename HType> class MultiHolder<1,HType> : public HType {
public:

    // allow implicit construction of MultiHolder<1,X> from X
    MultiHolder<1,HType>() : HType() {}
    MultiHolder<1,HType>(const HType& thing) : HType(thing) {}

    inline HType& operator[](size_t idx) { return *this; }
    inline HType& operator()(size_t idx) { return *this; }
    inline const HType& operator[](size_t idx) const { return *this; }
    inline const HType& operator()(size_t idx) const { return *this; }

};

// stream operator for BasePlaneWaveHolder class
template<int NUM_LAMBDAS, int NUM_WAVES>
std::ostream& operator<<(std::ostream& os, const BasePlaneWaveHolder<NUM_LAMBDAS, NUM_WAVES>& pw)
{
    os << pw.str();
    return os;
}

template<typename HType, int multiplicity>
inline std::ostream& operator<<(std::ostream& os, const MultiHolder<multiplicity,HType>& holder)
{
    os << holder.str();
    return os;
}

template<int NUM_TERMS>
class LegendreHolder : public TriangularMemory<NUM_TERMS, double> {};


#endif /* MULTIPOLE_HOLDER_H_ */
