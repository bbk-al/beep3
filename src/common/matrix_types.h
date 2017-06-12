/*
 * matrix_types.h
 *
 *  Created on: 16 Aug 2010
 *      Author: david
 */

#ifndef MATRIX_TYPES_H_
#define MATRIX_TYPES_H_
//#ifndef LTL_RANGE_CHECKING
//#define LTL_RANGE_CHECKING
//#endif
//#ifndef LTL_USE_SIMD
//#define LTL_USE_SIMD
//#endif

// These #defs cause annoying warnings when they
// collide with Charm
#ifdef __CHARMC__
#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif
#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif
#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif
#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif
#endif

#include <ltl/marray.h>
#include <ltl/fmatrix.h>
#include <ltl/fvector.h>

#include <complex>
#include <memory.h>

using ltl::MArray;
using ltl::Range;
using ltl::FMatrix;
using ltl::FVector;
using ltl::Shape;

typedef MArray<double, 1> DblMatrix1D;
typedef MArray<double, 2> DblMatrix2D;
typedef MArray<double, 3> DblMatrix3D;

typedef MArray<std::complex<double>, 1> CmplxMatrix1D;
typedef MArray<std::complex<double>, 2> CmplxMatrix2D;
typedef MArray<std::complex<double>, 3> CmplxMatrix3D;
typedef MArray<std::complex<double>, 4> CmplxMatrix4D;

typedef MArray<int, 1> UShrtMatrix1D;
typedef MArray<int, 2> UShrtMatrix2D;
typedef MArray<int, 3> UShrtMatrix3D;

// OK so this should really be an FMatrix to take advantage of the 
// LTL templatey goodness.  However when it gets serialised into
// Charm++ objects the FMatrix stuff goes belly-up (that is to say,
// pear-shaped) so I'm implementing it as a simple array which maps
// into a simple block of data in memory. (Because it has only the 
// array of 9 doubles, no virtual functions or anything else)
template <typename PrecisionType>
class Matrix_3x3
{
public:
    
    Matrix_3x3() {
        memset(data,0,sizeof(PrecisionType)*9);
    }
    Matrix_3x3(PrecisionType xx){
        for (int ii=0; ii < 9; ++ii)
        {
            data[ii] = xx;
        }
    }
    
    // zero-based indexing
    inline PrecisionType operator() (const int i, const int j) const
    {
        return data[i*3 + j];
    }
    
    // zero-based indexing
    inline PrecisionType& operator() (const int i, const int j)
    {
        return data[i*3 + j];
    
    }
    
    inline Matrix_3x3& operator+=(const Matrix_3x3& other)
    {
        for (int ii=0; ii < 9; ++ii)
        {
            data[ii] += other.data[ii];
        }   
        return *this;
    }
    
    template <typename OtherType>
    inline Matrix_3x3& operator+=(const OtherType* other)
    {
        for (int ii=0; ii < 9; ++ii)
        {
            data[ii] += static_cast<PrecisionType>(other[ii]);
        }   
        return *this;
    }
    
    inline void reset()
    {
        for (int ii=0; ii < 9; ++ii)
        {
            data[ii] = 0;
        }   
        return;
    }

protected:
    
    PrecisionType data[9];
    
};

typedef Matrix_3x3<double> GradField3x3;

#endif /* MATRIX_TYPES_H_ */
