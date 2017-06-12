/*
* fmm_time_info.h
*
*  Created on: 21 Jul 2010
*      Author: david
*/

#ifndef FMM_TIME_INFO_H_
#define FMM_TIME_INFO_H_

// you'll probably want to use this...
#include "../common/useful_clock.h"
#include <sstream>
#include <iostream>

namespace fmm
{

// This very simple class (could be a struct, but want to init. it zero, so need ctor)
// holds timings (return by the myclock() function) which can be incremented by 
// fMM functions to give very simple (and probably not massively precise) profile
// of timings.

class TimeInfo
{
private:
    
    // number of variables in the timing info
    static const int SIZE=13;
    
public:
    
    TimeInfo() 
    {
        
        // create array of pointers to data
        init_ptrs();
        
        // set all vars to zero
        zero();
    }
    
    TimeInfo(const TimeInfo& other) 
    {
        
        // create pointer array
        init_ptrs();
        
        // set vars to values contained in other
        for( int ii=0; ii < SIZE; ++ii)
        {
            *(this->data[ii]) = *(other.data[ii]);
        }
     }

    inline void init_ptrs() 
    {
        // create internal array of pointers to the actual data
        data[0] = &init;
        data[1] = &allocate_memory;
        data[2] = &create_mpoles;
        data[3] = &pass_mpoles_upwards;
        data[4] = &total_upward_pass;
        data[5] = &inherit_locals_downwards;
        data[6] = &interaction_list_xlations;
        data[7] = &convert_to_planewaves;
        data[8] = &convert_from_planewaves;
        data[9] = &total_downward_pass;
        data[10] = &fmm_evaluations;
        data[11] = &explicit_evaluations;
        data[12] = &num_ilist_xlations;
    }

    inline void zero() {
        
        // iterate over all items and set to zero
        for( int ii=0; ii < SIZE; ++ii)
        {
            *(data[ii]) = 0;
        }
    }

    inline TimeInfo& operator=(const TimeInfo& other)
    {
        // assign
        for( int ii=0; ii < SIZE; ++ii)
        {
            *(this->data[ii]) = *(other.data[ii]);
        }
        return *this;
    }

    inline TimeInfo operator+(const TimeInfo& other) const
    {
        TimeInfo new_blk(*this);
        new_blk += other;
        return new_blk;
    }

    inline TimeInfo operator/(int div) const
    {
        TimeInfo new_blk(*this);
        new_blk /= div;
        return new_blk;
    }

    inline TimeInfo operator*(int mult) const
    {
        TimeInfo new_blk(*this);
        new_blk *= mult;
        return new_blk;
    }

    inline TimeInfo& operator+=(const TimeInfo& other)
    {
        // add items
        for( int ii=0; ii < SIZE; ++ii)
        {
            *(this->data[ii]) += *(other.data[ii]);
        }
        return *this;
    }
    
    inline TimeInfo& operator*=(int mult)
    {
        // add items
        for( int ii=0; ii < SIZE; ++ii)
        {
            *(this->data[ii]) *= mult;
        }
        return *this;
    }
    
    inline TimeInfo& operator/=(int div)
    {
        // add items
        for( int ii=0; ii < SIZE; ++ii)
        {
            *(this->data[ii]) /= div;
        }
        return *this;
    }

    long init;
    long allocate_memory;
    
    long create_mpoles;
    long pass_mpoles_upwards;
    long total_upward_pass;
    
    long inherit_locals_downwards;
    long interaction_list_xlations;
    long convert_to_planewaves;
    long convert_from_planewaves;
    long total_downward_pass;
    
    long fmm_evaluations;
    long explicit_evaluations;
    long num_ilist_xlations;

    std::string str() const {
        
        std::ostringstream buf;
        for (int ii=0; ii < SIZE; ++ii)
        {
            if (ii > 0) { buf << ","; }
            buf << *(data[ii]);
        }
        return buf.str();   
    }
    
    static std::string describe() {
        
        std::ostringstream buf;
        buf << "init, "
            << "allocate_memory, "
            << "create_mpoles, "
            << "pass_mpoles_upwards, "
            << "total_upward_pass, "
            << "inherit_locals_downwards, "
            << "interaction_list_xlations, "
            << "convert_to_pw, "
            << "convert_from_pw, "
            << "total_downward_pass, "
            << "fmm_evals, "
            << "explicit_evals, "
            << "num_ilist_xlations\n";
        return buf.str();   
    }    
    
private:
    
    long* data[SIZE];

};

inline std::ostream& operator<<(std::ostream& os, const TimeInfo& t)
{
    os << t.str();
    return os;
}


} // end namespace


#endif // end #ifdef header guard



