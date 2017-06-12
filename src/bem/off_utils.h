/*
* off_utils.h
*
*  Created on: 21 Jul 2010
*      Author: david
*/
#ifndef _BEEP_OFF_UTILS_H_
#define _BEEP_OFF_UTILS_H_

#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include "triangle.h"
#include "vertex.h"
#include <exception>

// utility class for handy OFF-parsing functions
class OFF_Surface
{

public:

    static void read_mesh_from_file(const std::string& off_filename,
                                    std::vector<Vertex>& vertex_list,
                                    std::vector<Triangle>& triangle_list)
    {
        size_t line_ctr=0;
        
        // open the file
        std::ifstream off_file;
        off_file.open(off_filename.c_str());
        assert(off_file.good());

        // read in off-format data
        unsigned int num_verts, num_faces, num_edges;

        // skip comments
        while (off_file.peek() == '#' || off_file.peek() == '!') { ++line_ctr; off_file.ignore( std::numeric_limits<std::streamsize>::max(), '\n' ); assert(!off_file.eof()); }

        std::string first_line;
        getline(off_file, first_line);
        ++line_ctr;
        if (first_line != "OFF" && first_line != "nOFF") {
            
            //std::cerr << "Error reading OFF file " << off_filename << " -- first line should be OFF (optional flags not supported by BEEP)" << std::endl;
            try
            {
                std::stringstream buf;
                buf << first_line;
                buf >> num_verts >> num_faces >> num_edges;
                buf.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
            }
            catch (...)
            {
                std::cerr << "Error reading OFF file " << off_filename << " -- first line should read \"OFF\" or \"nOFF\"" << std::endl;
                std::cerr << "(Since OFF header keyword is missing, I expected \"num_verts num_faces num_edges\")" << std::endl;
                throw std::exception();
            }
            
        }
        else
        {
            off_file >> num_verts >> num_faces >> num_edges;
            off_file.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
            ++line_ctr;
        }

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

            // ignore comments
            while (off_file.peek() == '#') { 
                ++line_ctr; 
                off_file.ignore( std::numeric_limits<std::streamsize>::max(), '\n' ); 
                if(off_file.eof()) {
                    std::cerr << "Unexpected end of file " << off_filename << " at line " << line_ctr << ": Expected " << num_verts << " vertices, only found " << vertex_list.size() << std::endl;
                    throw std::exception();
                } 
            }
            
            // Could do with fixing this code repetition...
            if(off_file.eof()) {
                std::cerr << "Unexpected end of file " << off_filename << " at line " << line_ctr << ": Expected " << num_verts << " vertices, only found " << vertex_list.size() << std::endl;
                throw std::exception();
            } 

            std::string line;
            getline(off_file, line);
            ++line_ctr;
            std::stringstream buf(line);
            buf >> vx >> vy >> vz >> nx >> ny >> nz;
            if (nx != 0 || ny != 0 || nz != 0)
            {
                vnormals_set = true;
            }
            vertex_list.push_back(Vertex(Vector(vx,vy,vz), Vector(nx,ny,nz)));

        }
        //std::cout << "Added " << vertex_list.size() << " vertices." << std::endl;

        // get triangles
        unsigned int v1_idx, v2_idx, v3_idx, num;
        for (unsigned int tidx=0; tidx < num_faces; ++tidx)
        {
            while (off_file.peek() == '#') { 
                ++line_ctr; 
                off_file.ignore( std::numeric_limits<std::streamsize>::max(), '\n' ); 
                if(off_file.eof()) {
                    std::cerr << "Unexpected end of file " << off_filename << " at line " << line_ctr << ": Expected " << num_faces << " triangles, only found " << triangle_list.size() << std::endl;
                    throw std::exception();
                } 
            }
            
            if(off_file.eof()) {
                std::cerr << "Unexpected end of file " << off_filename << " at line " << line_ctr << ": Expected " << num_faces << " triangles, only found " << triangle_list.size() << std::endl;
                throw std::exception();
            } 
            off_file >> num >> v1_idx >> v2_idx >> v3_idx;
            ++line_ctr;
            if (num != 3) {
                std::cerr << "BEEP requires triangular meshes, so expected number of vertices per face to be 3 (" << off_filename << ":" << line_ctr << ")" << std::endl;
                throw std::exception();
            }
            off_file.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

            unsigned int t_idx = triangle_list.size();
            triangle_list.push_back(Triangle(vertex_list, v1_idx, v2_idx, v3_idx));

            vertex_list[v1_idx].add_triangle(t_idx);
            vertex_list[v2_idx].add_triangle(t_idx);
            vertex_list[v3_idx].add_triangle(t_idx);

        }
        //std::cout << "Added " << triangle_list.size() << " triangles." << std::endl;
        
        for (unsigned int vctr=0; vctr < vertex_list.size(); ++vctr)
        {
            Vertex& v = vertex_list[vctr];

            // reorder vertices to ensure contiguous triangles
            v.reorder(triangle_list);
            
            const std::vector<unsigned int>& triangle_indices = v.get_triangle_indices();
            
            for (std::vector<unsigned int>::const_iterator it=triangle_indices.begin(), next, end=triangle_indices.end(); it != end; ++it)
            {
                next = it+1;
                if (next == end) { next = triangle_indices.begin(); }
                
                Triangle& t1 = triangle_list[*it];
                Triangle& t2 = triangle_list[*next];
                t1.set_adjacent_normal(t2);
                t2.set_adjacent_normal(t1);
            }
            
            v.set_normal(triangle_list);
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
            //std::cout << "Vertex normals were read in with the off file." << std::endl;
        }

        return;
    }
};

#endif // _BEEP_OFF_UTILS_H_
