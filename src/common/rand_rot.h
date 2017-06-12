#ifndef __RAND_ROT_H__
#define __RAND_ROT_H__

#include "math_vector.h" // for Quaternion def
#include <time> // for seeding the rng
#include <boost/random/uniform_real.hpp> // for random numbers in range zero to 1
#include <math.h> // for sqrt

// Random rotation- is a random point on a 4-dimensional unit sphere.
// See Numerical Recipes. (Press et al.)
class RandomRotation : public Quaternion
{
    RandomRotation() : Quaternion(0,0,0,0)
    {
        generator.seed(static_cast<unsigned int>(std::time(0)));    
        boost::uniform_real<> uni_dist(-1,1);
        boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);
        
        double u0=0,u1=0,u2=0,u3=0;
        do
        {
            u0 = uni();
            u1 = uni();
        } while (u0*u0 + u1*u1 > 1.0)
        do
        {
            u2 = uni();
            u3 = uni();
        } while (u2*u2 + u3*u3 > 1.0)
        
        this->a = u0;
        this->b = u1;
        double useful = sqrt((1.0 - u0*u0 - u1*u1)/(u2*u2 + u3*u3));
        this->c = u2 * useful;
        this->d = u3 * useful;
        
        return;

    }
    
};



#endif // __RAND_ROT_H__
