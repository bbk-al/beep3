/*
 * octree_indexer.h
 *
 *  Created on: 16 Aug 2010
 *      Author: david
 */

#ifndef OCTREE_INDEXER_H_
#define OCTREE_INDEXER_H_

#include <cassert>

#include "math_vector.h"
#include "matrix_types.h"
#include "multipole_holder.h"

#ifdef __INTEL_COMPILER
#include <mathimf.h>
#else
#include <math.h>
#endif
#include <sstream>
#include <iostream>
#include <string>
#include <complex>
#include <limits>
#include <memory>
#include <algorithm>
#include <vector>

#include <iterator>
#include <memory>
#include <boost/shared_array.hpp>

#define PRECIS std::numeric_limits<double>::epsilon()
#define MIN_HACKY(x,y) ((x<y) ? x : y)

static const unsigned long two_powers[64] = {1ul<<0, 1ul<<1, 1ul<<2, 1ul<<3, 1ul<<4, 1ul<<5, 1ul<<6, 1ul<<7,
		                                        1ul<<8, 1ul<<9, 1ul<<10, 1ul<<11, 1ul<<12, 1ul<<13, 1ul<<14, 1ul<<15,
		                                        1ul<<16, 1ul<<17, 1ul<<18, 1ul<<19, 1ul<<20, 1ul<<21, 1ul<<22, 1ul<<23,
		                                        1ul<<24, 1ul<<25, 1ul<<26, 1ul<<27, 1ul<<28, 1ul<<29, 1ul<<30, 1ul<<31,
		                                        1ul<<32, 1ul<<33, 1ul<<34, 1ul<<35, 1ul<<36, 1ul<<37, 1ul<<38, 1ul<<39,
		                                        1ul<<40, 1ul<<41, 1ul<<42, 1ul<<43, 1ul<<44, 1ul<<45, 1ul<<46, 1ul<<47,
		                                        1ul<<48, 1ul<<49, 1ul<<50, 1ul<<51, 1ul<<52, 1ul<<53, 1ul<<54, 1ul<<55,
		                                        1ul<<56, 1ul<<57, 1ul<<58, 1ul<<59, 1ul<<60, 1ul<<61, 1ul<<62, 1ul<<63};

class BadIndexer : public std::exception
{
    public:
        BadIndexer() : std::exception() {};
};

class OctreeIndexer
{

    // This class is to create a unique hash key for Octree nodes --
    // it must incorporate the level in the tree to which the node
    // belongs, plus the x/y/z index of the node, all within the size
    // of 3 ints.

public:

	static const unsigned short MAX_LEVELS = MIN_HACKY((sizeof(int)*8) - 2, (sizeof(unsigned long)*8 / 3)-1);
    unsigned int data[3];

    OctreeIndexer()
    {
        data[0]=0;
        data[1]=0;
        data[2]=0;
    };

    OctreeIndexer(const OctreeIndexer& other)
    {
        //memcpy(data, other.data, sizeof(data));
        data[0] = other.data[0];
        data[1] = other.data[1];
        data[2] = other.data[2];
        return;
    };

    OctreeIndexer(unsigned int level, unsigned int x_idx, unsigned int y_idx, unsigned int z_idx)
    {
        init(level, x_idx, y_idx, z_idx, data);

        // check we can unpack what we think we packed...
        assert(get_x_idx() == x_idx);
        assert(get_y_idx() == y_idx);
        assert(get_z_idx() == z_idx);
        assert(get_level() == level);

        return;
    }

    inline unsigned short id_within_parent() const
    {
        unsigned int x_idx = get_x_idx();
        unsigned int y_idx = get_y_idx();
        unsigned int z_idx = get_z_idx();

        unsigned short id = 1;
        if (x_idx % 2) {id += 1;}
        if (y_idx % 2) {id += 2;}
        if (z_idx % 2) {id += 4;}
        return id;
    }

    inline bool operator<(const OctreeIndexer& other) const
    {
        return this->as_hash_number() < other.as_hash_number();
    }

    inline unsigned long two_pow(const unsigned char power) const
    {
    	return two_powers[power];
    }

    inline unsigned long __hash() const
    {
        unsigned short level = get_level();
        
        unsigned long hash = 0;
        for (unsigned short ii=0; ii < level; ++ii)
        {
            hash += two_pow(ii*3);
        }
        
        OctreeIndexer xx = *this;
        unsigned long mult=1;
        while (level)
        {
            hash += (xx.id_within_parent()-1)*mult;
            mult *= 8;
            xx = xx.translate_upwards_to_level(level-1);
            --level;
        }  
        
        return hash;
        
    }

    inline unsigned long as_hash_number() const
    {
        // the hash number gives the unique id of this cell within the mxmxm cube
        // where m=2^n (n is the level number, which is the same thing as the
        // number of subdivisions of the root level cube)
#if 0
        // TODO: assert that an unsigned long is big enough
        return __hash();
#else
        unsigned short level = get_level();

        unsigned int x = get_x_idx();
        unsigned int y = get_y_idx();
        unsigned int z = get_z_idx();

        unsigned long int edge_cells = two_pow(level);

        unsigned long hash = 0;
        for (unsigned short ii=0; ii < level; ++ii)
        {
            hash += two_pow(ii*3);
        }
        /*
        // get sixteenmer block
        unsigned long int edge_cells = (level >= 4) ? pow(2,level-4) : 0;
        hash += ((x / 16) + (y / 16)*edge_cells + (z / 16)*edge_cells*edge_cells)*16*16*16;

        // get remainders
        x %= 16;
        y %= 16;
        z %= 16;

        // get index of octomer
        hash += (x/2 + (y/2)*8 + (z/2)*64)*8;

        // get index within octet
        hash += id_within_parent();*/

        hash += x;
        hash += y*edge_cells;
        hash += z*edge_cells*edge_cells;

        return hash;
#endif
    }

    static void init(unsigned int level, unsigned int x_idx, unsigned int y_idx, unsigned int z_idx, unsigned int data[3])
    {
        // Highly shonky!!  I'm using the three int fields for x,y,z indices, but left-shifting each
        // by two bits (i.e. multiply by 4) to make room to store the level in the 6 bits this creates.
        // Therefore the level must be <= 63 (maximum unsigned number with 6 bits), and the x/y/z indices
        // must have at least 2 bits spare at their upper end -- i.e. be less than a quarter of the maximum
        // size of an unsigned int.
        // I'm using multiplication rather than explicit bit shifting to avoid endian-ness problems and
        // non-portability.  (Wear your extra tight corset whilst reading this).

        // cannot handle level being more than 6 bits
        assert(level <= 30);

//         data[0]=0;
//         data[1]=0;
//         data[2]=0;

        // create 2 bits of room at the least-significant end
        unsigned int shifted_x_idx = x_idx * 4;
        unsigned int shifted_y_idx = y_idx * 4;
        unsigned int shifted_z_idx = z_idx * 4;

        // check that we didn't just overflow those variables
        if(    !(((shifted_x_idx / 4) == x_idx)
                && ((shifted_y_idx / 4) == y_idx)
                && ((shifted_z_idx / 4) == z_idx)) )
		{
			std::cout << x_idx << " " << y_idx << " " << z_idx << "\n";
			std::cout << shifted_x_idx << " " << shifted_y_idx << " " << shifted_z_idx << "\n";
		}


        assert(    ((shifted_x_idx / 4) == x_idx)
                && ((shifted_y_idx / 4) == y_idx)
                && ((shifted_z_idx / 4) == z_idx));

        // a bit of bit fiddling here
        if (level & 0x01) { shifted_x_idx += 1; }
        if (level & 0x02) { shifted_x_idx += 2; }
        if (level & 0x04) { shifted_y_idx += 1; }
        if (level & 0x08) { shifted_y_idx += 2; }
        if (level & 0x10) { shifted_z_idx += 1; }
        if (level & 0x20) { shifted_z_idx += 2; }

        data[0] = shifted_x_idx;
        data[1] = shifted_y_idx;
        data[2] = shifted_z_idx;

        return;

    }

    inline bool operator==(const OctreeIndexer& other) const {
        return (   data[0] == other.data[0]
                && data[1] == other.data[1]
                && data[2] == other.data[2]);
    }

    inline bool operator!=(const OctreeIndexer& other) const {
        return (   data[0] != other.data[0]
                || data[1] != other.data[1]
                || data[2] != other.data[2]);
    }

    inline unsigned int get_x_idx() const
    {
        return data[0] / 4;
    }

    inline unsigned int get_y_idx() const
    {
        return data[1] / 4;
    }

    inline unsigned int get_z_idx() const
    {
        return data[2] / 4;
    }

    inline unsigned short get_level() const
    {
        unsigned short level=0;
        if (data[0] & 0x01) { level += 1; }
        if (data[0] & 0x02) { level += 2; }
        if (data[1] & 0x01) { level += 4; }
        if (data[1] & 0x02) { level += 8; }
        if (data[2] & 0x01) { level += 16; }
        if (data[2] & 0x02) { level += 32; }
        return level;
    }

    // function to return an OctreeIndexer object
    // which refers to the child Octant of that referred
    // to by this OctreeIndexer.
    inline OctreeIndexer get_child_idx(int child_id) const
    {
        assert(child_id >= 1 && child_id <= 8);

        unsigned short child_level = get_level() + 1;

        unsigned int child_x_idx = get_x_idx() * 2;
        unsigned int child_y_idx = get_y_idx() * 2;
        unsigned int child_z_idx = get_z_idx() * 2;

        if (child_id % 2 == 0) { child_x_idx++; }
        if ((child_id % 4) == 0 || (child_id % 4) == 3) { child_y_idx++; }
        if (child_id > 4) { child_z_idx++; }

        return OctreeIndexer(child_level, child_x_idx, child_y_idx, child_z_idx);
    }

    inline OctreeIndexer get_parent_idx() const
    {
        return OctreeIndexer(get_level() - 1, get_x_idx() / 2, get_y_idx() / 2, get_z_idx() / 2);
    }

    inline OctreeIndexer translate_upwards_to_level(unsigned short target_level) const
    {
        unsigned short cur_level = get_level();
        assert(cur_level >= target_level);
        unsigned int div = static_cast<unsigned int>(two_pow(cur_level - target_level));
        return OctreeIndexer(target_level, get_x_idx() / div, get_y_idx() / div, get_z_idx() / div);
    }

    inline bool is_root() const { return (get_level() == 0); }

    inline bool is_up()    const { return (get_z_idx() % 2); }
    inline bool is_north() const { return (get_y_idx() % 2); }
    inline bool is_east()  const { return (get_x_idx() % 2); }

    inline void get_transmission_neighbours(unsigned short chare_level, std::vector<OctreeIndexer>& transmit_here) const
    {

        // figure out what chare this is
        OctreeIndexer this_chare = this->translate_upwards_to_level(chare_level);

        for (short nb_id=0; nb_id < 27; ++nb_id)
        {
            try
            {
                OctreeIndexer neighbour(get_neighbour_idxer(nb_id));
                OctreeIndexer nb_chare(neighbour.translate_upwards_to_level(chare_level));
                if (nb_chare != this_chare)
                {
                    // add neighbour chare to the list (if it's not already there)
                    std::vector<OctreeIndexer>::const_iterator find_it = std::find(transmit_here.begin(), transmit_here.end(), nb_chare);
                    if (find_it == transmit_here.end()) {
                        transmit_here.push_back(nb_chare);
                    }
                }
            }
            catch (BadIndexer) {}
        }
    }

    std::string str() const
    {
        std::ostringstream oss;
        oss << "( level: " << get_level() << " x/y/z: (" << get_x_idx() << "," << get_y_idx() << "," << get_z_idx() << ") [hash=" << as_hash_number() << "]";
        return oss.str();
    }

    OctreeIndexer get_neighbour_idxer(unsigned short neighbour_id) const
    {
        assert(neighbour_id <= 26); // max number of neighbours is 27
        if (neighbour_id == 13) { return *this; }

        int x_idx = static_cast<int>(get_x_idx());
        int y_idx = static_cast<int>(get_y_idx());
        int z_idx = static_cast<int>(get_z_idx());
        const unsigned short level = get_level();
        const int edge_cells_at_level = static_cast<int>(two_pow(level));

        if (neighbour_id >= 18) {
            z_idx++;
        } else if (neighbour_id <= 8) {
            z_idx--;
        }
        if (z_idx < 0 || z_idx >= edge_cells_at_level) { throw BadIndexer(); }

        neighbour_id %= 9;

        if (neighbour_id >= 6) {
            y_idx++;
        } else if (neighbour_id <= 2) {
            y_idx--;
        }
        if (y_idx < 0 || y_idx >= edge_cells_at_level) { throw BadIndexer(); }

        neighbour_id %= 3;

        if (neighbour_id == 0) {
            x_idx--;
        }
        else if (neighbour_id == 2) {
            x_idx++;
        }
        if (x_idx < 0 || x_idx >= edge_cells_at_level) { throw BadIndexer(); }

        return OctreeIndexer(level, x_idx, y_idx, z_idx);
    }

};

// stream operator for CharmNodePatchIndex class
inline std::ostream& operator<<(std::ostream& os, const OctreeIndexer& idxer)
{
    os << idxer.str();
    return os;
}

#endif /* OCTREE_INDEXER_H_ */
