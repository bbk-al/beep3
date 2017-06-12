/*
* gts_utils.h
*
*  Created on: 21 Jul 2010
*      Author: david
*/
#ifndef _BEEP_GTS_UTILS_H_
#define _BEEP_GTS_UTILS_H_

#include <algorithm>
#include <vector>
#include <set>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <istream>
#include <limits>
#include "triangle.h"
#include "vertex.h"
#include "edge.h"

// utility class for handy GTS-parsing functions
class GTS_Surface
{

public:

    static void read_mesh_from_file(const std::string& gts_filename,
                                    std::vector<Vertex>& vertex_list,
                                    std::vector<Triangle>& triangle_list)
    {
        std::vector<Edge> edge_list;
        
        // open the file
        std::ifstream gts_file;
        gts_file.open(gts_filename.c_str());
        assert(gts_file.good());

        // read in gts-format data
        unsigned int num_verts, num_edge_list, num_faces;

        // skip comments
        while (gts_file.peek() == '#' || gts_file.peek() == '!') { gts_file.ignore( std::numeric_limits<std::streamsize>::max(), '\n' ); assert(!gts_file.eof()); }

        gts_file >> num_verts >> num_edge_list >> num_faces;
        gts_file.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
        assert(num_edge_list == (3*num_faces / 2));

        // reserve storage space
        typedef std::set<unsigned int> EdgeSet;
        EdgeSet* gts_edge_list = new EdgeSet[num_edge_list];
        assert(gts_edge_list != NULL);
        edge_list.reserve(num_edge_list);
        vertex_list.reserve(num_verts);
        triangle_list.reserve(num_faces);

        // get vertices
        bool vnormals_set = false;
        double vx,vy,vz;
        double nx=0,ny=0,nz=0;
        for (unsigned int vidx=0; vidx < num_verts; ++vidx)
        {
            nx=0;
            ny=0;
            nz=0;

            while (gts_file.peek() == '#' || gts_file.peek() == '!') { gts_file.ignore( std::numeric_limits<std::streamsize>::max(), '\n' ); assert(!gts_file.eof()); }
            std::string line;
            getline(gts_file, line);
            std::stringstream buf(line);
            buf >> vx >> vy >> vz >> nx >> ny >> nz;
            if (nx != 0 || ny != 0 || nz != 0)
            {
                vnormals_set = true;
            }
            vertex_list.push_back(Vertex(Vector(vx,vy,vz), Vector(nx,ny,nz)));

        }
        //std::cout << "Added " << vertex_list.size() << " vertices." << std::endl;

        // get edge_list
        unsigned int edge_vert_1, edge_vert_2;
        for (unsigned int eidx=0; eidx < num_edge_list; ++eidx)
        {
            while (gts_file.peek() == '#' || gts_file.peek() == '!') { gts_file.ignore( std::numeric_limits<std::streamsize>::max(), '\n' ); assert(!gts_file.eof()); }
            gts_file >> edge_vert_1 >> edge_vert_2;
            gts_file.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
            gts_edge_list[eidx].insert(edge_vert_1 - 1);
            gts_edge_list[eidx].insert(edge_vert_2 - 1);
            edge_list.push_back(Edge());
            edge_list[eidx].set_vertices(edge_vert_1 - 1, edge_vert_2 - 1);
        }

        // get triangles
        unsigned int a,b,c;
        for (unsigned int tidx=0; tidx < num_faces; ++tidx)
        {
            while (gts_file.peek() == '#' || gts_file.peek() == '!') { gts_file.ignore( std::numeric_limits<std::streamsize>::max(), '\n' ); assert(!gts_file.eof()); }
            gts_file >> a >> b >> c;
            gts_file.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

            std::set<unsigned int> intersection;
            std::set<unsigned int>& e1 = gts_edge_list[a-1];
            std::set<unsigned int>& e2 = gts_edge_list[b-1];
            std::set<unsigned int>& e3 = gts_edge_list[c-1];

            intersection.clear();
            std::set_intersection(e1.begin(),e1.end(),e2.begin(),e2.end(), std::inserter(intersection, intersection.begin()));
            unsigned int v1_idx = *(intersection.begin());

            intersection.clear();
            std::set_intersection(e2.begin(),e2.end(),e3.begin(),e3.end(), std::inserter(intersection, intersection.begin()));
            unsigned int v2_idx = *(intersection.begin());

            intersection.clear();
            std::set_intersection(e3.begin(),e3.end(),e1.begin(),e1.end(), std::inserter(intersection, intersection.begin()));
            unsigned int v3_idx = *(intersection.begin());

            unsigned int t_idx = triangle_list.size();
            triangle_list.push_back(Triangle(vertex_list, v1_idx,v2_idx,v3_idx));

            vertex_list[v1_idx].add_triangle(t_idx);
            vertex_list[v2_idx].add_triangle(t_idx);
            vertex_list[v3_idx].add_triangle(t_idx);


            edge_list[a-1].add_triangle(tidx);
            edge_list[b-1].add_triangle(tidx);
            edge_list[c-1].add_triangle(tidx);
        }
        //std::cout << "Added " << triangle_list.size() << " triangles." << std::endl;

        for (std::vector<Edge>::const_iterator it=edge_list.begin(), end=edge_list.end(); it != end; ++it)
        {
            // assert that all edge_list have exactly two triangles (otherwise mesh is pathological)
            assert(it->tctr == 2);
            const Edge& ed = *it;
            Triangle& t1 = triangle_list[ed.t1_idx];
            Triangle& t2 = triangle_list[ed.t2_idx];
            t1.set_adjacent_normal(t2);
            t2.set_adjacent_normal(t1);
        }

        if (vnormals_set == false)
        {
            //std::cout << "Calculating vertex normals by mean planar weight around nodes." << std::endl;
            for (std::vector<Vertex>::iterator v_it=vertex_list.begin(), v_end=vertex_list.end(); v_it != v_end; ++v_it)
            {
                v_it->set_normal(triangle_list);
            }

            // reset the vertex normals within the triangle object
            for (std::vector<Triangle>::iterator tri_it=triangle_list.begin(), tri_end=triangle_list.end(); tri_it != tri_end; ++tri_it)
            {
                tri_it->reinit_vertex_normals(vertex_list);
            }
        }
        else
        {
            // verify that all vertex normals are set and valid
            for (std::vector<Vertex>::const_iterator vit=vertex_list.begin(), vend=vertex_list.end(); vit != vend; ++vit)
            {
                const Vertex& v = *vit;
                assert(fabs(v.get_normal().length() - 1.0) < 1e-6);
            }
            //std::cout << "Vertex normals were read in with the gts file." << std::endl;
        }

        // free allocated memory
        delete[] gts_edge_list;

        return;
    }
};

#endif // _BEEP_GTS_UTILS_H_
