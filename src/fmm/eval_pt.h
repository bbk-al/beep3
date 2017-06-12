/*
* eval_pt.h
*
*  Created on: 22 Sep 2010
*      Author: david
*/

#ifndef EVAL_PT_H_
#define EVAL_PT_H_

#include <boost/utility.hpp>
#define THREADSAFE // comment this out if not paranoid about thread safety

// If using OpenCL, MUST use threadsafe evaluation points
#ifdef OPENCL
#ifndef THREADSAFE
#define THREADSAFE
#endif
#endif

#ifdef THREADSAFE
#include <boost/thread/mutex.hpp>
#endif

#include "../common/math_vector.h"
#include "../common/matrix_types.h"
#include "../common/charge.h"

namespace fmm
{

// Handy little struct to represent an evaluation point in the FMM
class EvalPtBase : public Vector
{
public:

#ifdef THREADSAFE   
    boost::mutex mutex;
#endif
    
public:

    EvalPtBase() : Vector(0,0,0), unique_index(0), field(0,0,0), pot(0)  {}
    EvalPtBase(const EvalPtBase& other) : Vector(other), unique_index(other.unique_index), pot(other.pot), field(other.field) {}
    EvalPtBase(const Vector& where) : Vector(where), unique_index(0), pot(0), field(0,0,0) {}

    virtual ~EvalPtBase() {}

    virtual void init() {
#ifdef THREADSAFE   
        boost::mutex* in_place = new(&mutex) boost::mutex;
#endif
        Vector::reset();
        unique_index = 0;
        field.reset();
        pot = 0;
    }
    
    inline void reinit_mutex() {
#ifdef THREADSAFE   
        boost::mutex* in_place = new(&mutex) boost::mutex;
#endif
    }

    // combine eval points
    EvalPtBase& operator+=(EvalPtBase& other) {

        add_potential(other.pot);
        add_field(other.field);
        return *this;
    }

    virtual void add_raw(const float* results)
    {
#ifdef THREADSAFE
        boost::mutex::scoped_lock l(mutex);
#endif
        pot += results[0];
        field.x += results[1];
        field.y += results[2];
        field.z += results[3];
    }

    // assignment operator changes the location of the evaluation point
    EvalPtBase& operator=(const Vector& where) {
#ifdef THREADSAFE
        boost::mutex::scoped_lock l(mutex);
#endif
        static_cast<Vector&>(*this) = where;
        return *this;
    }
    
    // assignment operator (analagous to copy construct)
    EvalPtBase& operator=(const EvalPtBase& other) 
    {
#ifdef THREADSAFE
        boost::mutex::scoped_lock l(mutex);
#endif
        static_cast<Vector&>(*this) = static_cast<const Vector&>(other);
        unique_index = other.get_idx();
        pot = other.get_potential();
        field = other.get_field();
        return *this;
    }

    virtual void reset() {
#ifdef THREADSAFE       
        boost::mutex::scoped_lock l(mutex);
#endif
        pot = 0;
        field.x = 0;
        field.y = 0;
        field.z = 0;
    }

    // alias pt as the Vector facet of self
    inline Vector& pt() { return static_cast<Vector&>(*this); }
    inline const Vector& pt() const { return static_cast<const Vector&>(*this); }

    inline void set_idx(size_t idx) { unique_index = idx; }
    inline size_t get_idx() const { return unique_index; }

    // these get/set functions are controlled by mutex to ensure no
    // shonky thread-unsafe action
    inline Vector get_field() {
#ifdef THREADSAFE       
        boost::mutex::scoped_lock l(mutex);
#endif
        return field;
    }

    // non-thread safe version
    inline const Vector& get_field() const { return field; }

    inline void set_field(const Vector& val)
    {
#ifdef THREADSAFE       
        boost::mutex::scoped_lock l(mutex);
#endif
        field = val;

    }
    inline void add_field(const Vector& val)
    {
#ifdef THREADSAFE
        boost::mutex::scoped_lock l(mutex);
#endif
        field += val;
    }

    inline double get_potential() {
#ifdef THREADSAFE
        boost::mutex::scoped_lock l(mutex);
#endif
        return pot;
    }
    inline const double& get_potential() const { return pot; }

    inline void set_potential(double val)
    {
#ifdef THREADSAFE
        boost::mutex::scoped_lock l(mutex);
#endif
        pot = val;
    }
    inline void add_potential(double val)
    {
#ifdef THREADSAFE
        boost::mutex::scoped_lock l(mutex);
#endif
        pot += val;
    }

    // a public function to allow test program to use it
    static void add_explicit_contrib(const Vector& pt, const Charge& ch, double& pot, Vector& field, double beta)
    {
        const double rx = (pt.x - ch.x);
        const double ry = (pt.y - ch.y);
        const double rz = (pt.z - ch.z);
        const double d_squared = rx*rx + ry*ry + rz*rz;
        const double d = sqrt(d_squared);

        // don't bother if this (point) charge is at the evaluation position
        if (d_squared <= PRECIS){ return; }

        // add screened coulomb potential
        const double pot_term = ch.get_charge()*exp(-beta*d) / d;
        pot += pot_term;

        // first derivative of potential (field)
        const double field_term = - pot_term * (1.0 + d*beta) / d_squared;
        field.x += field_term * rx;
        field.y += field_term * ry;
        field.z += field_term * rz;

        return;
    }

    static void add_explicit_contrib(EvalPtBase& pt, const Charge& ch, double beta)
    {
        double pot=0;
        Vector field(0,0,0);
        add_explicit_contrib(pt, ch, pot, field, beta);
        pt.add_potential(pot);
        pt.add_field(field);
    }


protected:

    size_t unique_index;
    double pot;
    Vector field;
    

};

inline std::ostream& operator<<(std::ostream& os, const EvalPtBase& ep)
{
    os << "uid=" << ep.get_idx() << " @ (" << static_cast<const Vector&>(ep) << ")"
        << " pot="  << ep.get_potential()
        << " field=" << ep.get_field();

    return os;

}

class EvalPt_2ndDerivs : public EvalPtBase, public GradField3x3
{

public:

#ifdef THREADSAFE   
    boost::mutex mutex;
#endif
    
public:

    EvalPt_2ndDerivs() : EvalPtBase(), GradField3x3(0.0) { }
    EvalPt_2ndDerivs(const EvalPt_2ndDerivs& other) : EvalPtBase(other), GradField3x3(static_cast<const GradField3x3&>(other)) {}
    explicit EvalPt_2ndDerivs(const Vector& where) : EvalPtBase(where), GradField3x3(0.0) {}
    virtual ~EvalPt_2ndDerivs() {}    

    void init() {
#ifdef THREADSAFE   
        boost::mutex* in_place = new(&mutex) boost::mutex;
#endif
        EvalPtBase::init();
        unique_index = 0;
        GradField3x3::reset();
    }
    
    inline void reinit_mutex() {
#ifdef THREADSAFE   
        boost::mutex* in_place = new(&mutex) boost::mutex;
#endif
    }

    // combine eval points
    EvalPt_2ndDerivs& operator+=(EvalPt_2ndDerivs& other) {

        assert(unique_index == other.unique_index);
        add_potential(other.pot);
        add_field(other.field);
        add_field2(other);
        return *this;
    }

    // assignment operator changes the location of the evaluation point
    EvalPt_2ndDerivs& operator=(const Vector& where) {
#ifdef THREADSAFE
        boost::mutex::scoped_lock l(mutex);
#endif
        static_cast<Vector&>(*this) = where;
        return *this;
    }
    
    // assignment operator (analagous to copy construct)
    EvalPt_2ndDerivs& operator=(const EvalPt_2ndDerivs& other) 
    {
#ifdef THREADSAFE
        boost::mutex::scoped_lock l(mutex);
#endif
        static_cast<Vector&>(*this) = static_cast<const Vector&>(other);
        unique_index = other.get_idx();
        pot = other.get_potential();
        field = other.get_field();
        static_cast<GradField3x3&>(*this) = static_cast<const GradField3x3&>(other);
        return *this;
    }

    void reset() {
#ifdef THREADSAFE       
        boost::mutex::scoped_lock l(mutex);
#endif
        pot = 0;
        field.x = 0;
        field.y = 0;
        field.z = 0;
        GradField3x3::reset();
    }

    // these get/set functions are controlled by mutex to ensure no
    // shonky thread-unsafe action
    inline GradField3x3 get_field2() {
#ifdef THREADSAFE
        boost::mutex::scoped_lock l(mutex);
#endif
        return *this;
    }
    inline const GradField3x3& get_field2() const { return *this; }

    inline void set_field2(const GradField3x3& val)
    {
#ifdef THREADSAFE
        boost::mutex::scoped_lock l(mutex);
#endif
        static_cast<GradField3x3&>(*this) = val;

    }
    inline void add_field2(const GradField3x3& val)
    {
#ifdef THREADSAFE
        boost::mutex::scoped_lock l(mutex);
#endif
        static_cast<GradField3x3&>(*this) += val;
    }

    virtual void add_raw(const float* results)
    {
#ifdef THREADSAFE
        boost::mutex::scoped_lock l(mutex);
#endif
        pot += results[0];
        field.x += results[1];
        field.y += results[2];
        field.z += results[3];
        static_cast<GradField3x3&>(*this) += &(results[4]); // should add all 9 numbers from results[4] onwards
    }

    // a public function to allow test program to use it
    static void add_explicit_contrib(const Vector& pt, const Charge& ch, double& pot, Vector& field, GradField3x3& field2, double beta)
    {
        const double rx = (pt.x - ch.x);
        const double ry = (pt.y - ch.y);
        const double rz = (pt.z - ch.z);
        const double d_squared = rx*rx + ry*ry + rz*rz;
        const double d = sqrt(d_squared);

        // don't bother if this (point) charge is at the evaluation position
        if (d_squared <= PRECIS) {return;}

        // add screened coulomb potential
        const double pot_term = ch.get_charge()*exp(-beta*d) / d;
        pot += pot_term;

        // first derivative of potential (field)
        const double field_term = - pot_term * (1.0 + d*beta) / d_squared;
        field.x += field_term * rx;
        field.y += field_term * ry;
        field.z += field_term * rz;

        // diagonal terms
        const double dbeta = d*beta;
        const double d_beta_squared = dbeta*dbeta;
        const double d4 = d_squared * d_squared;
        field2(0,0) += -pot_term*( (1.0+dbeta)*(d_squared - 3.0*rx*rx) - rx*rx*d_beta_squared) / d4;
        field2(1,1) += -pot_term*( (1.0+dbeta)*(d_squared - 3.0*ry*ry) - ry*ry*d_beta_squared) / d4;
        field2(2,2) += -pot_term*( (1.0+dbeta)*(d_squared - 3.0*rz*rz) - rz*rz*d_beta_squared) / d4;

        // off-diagonals
        const double f2_mixed_term = pot_term*(3.0 + 3.0*dbeta + d_beta_squared)/d4;
        field2(0,1) += f2_mixed_term * rx*ry;
        field2(0,2) += f2_mixed_term * rx*rz;
        field2(1,0) += f2_mixed_term * ry*rx;
        field2(1,2) += f2_mixed_term * ry*rz;
        field2(2,0) += f2_mixed_term * rz*rx;
        field2(2,1) += f2_mixed_term * rz*ry;

        return;
    }

    static void add_explicit_contrib(EvalPt_2ndDerivs& pt, const Charge& ch, double beta)
    {
        double pot=0;
        Vector field(0,0,0);
        GradField3x3 field2;
        add_explicit_contrib(pt, ch, pot, field, field2, beta);
        pt.add_potential(pot);
        pt.add_field(field);
        pt.add_field2(field2);
        
    }

};

inline std::ostream& operator<<(std::ostream& os, const EvalPt_2ndDerivs& ep)
{
    GradField3x3 f2 = const_cast<EvalPt_2ndDerivs&>(ep).get_field2();

    os << "uid=" << ep.get_idx() << " @ (" << static_cast<const Vector&>(ep) << ")"
        << " pot="  << ep.get_potential()
        << " field=" << ep.get_field()
        << " field2={";
    for (int ii=0; ii < 3; ++ii)
    {
        os << "(" << f2(ii,0) << "," << f2(ii,1) << "," << f2(ii,2) << ")";
    }
    os << "}";

    return os;

}

class EvalPt : public EvalPtBase
{
    // Prevent confusion between EvalPt and EvalPt_2ndDeriv -- common base,
    // but one not directly derived from the other

public:

    EvalPt() : EvalPtBase() { }
    EvalPt(const EvalPt& other) : EvalPtBase(other) {}
    explicit EvalPt(const Vector& where) : EvalPtBase(where) {}
    virtual ~EvalPt() {}


};

} // end namespace

#endif /* EVAL_PT_H_ */
