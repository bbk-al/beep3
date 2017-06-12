/*
 * edge.h
 *
 *  Created on: 23 Jul 2010
 *      Author: david
 */

#ifndef EDGE_H_
#define EDGE_H_

#include <cassert>

class Edge
{
public:

    Edge() : v1_idx(0), v2_idx(0), t1_idx(0), t2_idx(0), tctr(0) {}
    virtual ~Edge() {}

    void set_vertices(unsigned int idx1, unsigned int idx2)
    {
        assert(idx1 != idx2);
        if (idx1 > idx2) {
            v1_idx = idx2;
            v2_idx = idx1;
        }
        else
        {
            v1_idx = idx1;
            v2_idx = idx2;
        }

        return;
    }
    void add_triangle(unsigned int idx)
    {
        assert(tctr < 2);

        if (tctr==0)
        {
            t1_idx = idx;
        }
        else
        {
            if (t1_idx > idx)
            {
                t2_idx = t1_idx;
                t1_idx = idx;
            }
            else
            {
                t2_idx = idx;
            }
        }

        tctr++;
    }
#ifdef __CHARMC__
    void pup(PUP::er &p) {

        p | v1_idx;
        p | v2_idx;
        p | t1_idx;
        p | t2_idx;
    }
#endif

    unsigned int v1_idx;
    unsigned int v2_idx;
    unsigned int t1_idx;
    unsigned int t2_idx;
    unsigned short tctr;
};

#endif /* EDGE_H_ */
