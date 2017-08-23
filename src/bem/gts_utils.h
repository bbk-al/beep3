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
#include <memory>
#include "vertex.h"

// utility class for handy GTS-parsing functions
class GTS_Surface
{

public:

	template <unsigned int N>
	using IdxListN = std::vector<std::array<unsigned int,N>>;
	// (shared_ptr to match Mesh definition, unique would be ok)

    static void read_mesh_from_file(
		const std::string& gts_filename,	// in, file to read
		std::vector<Vertex>& vertex_list,	// out, vertices with normals set
		IdxListN<2>& gts_edge_list,		// out, vertex index pairs, low-high
		IdxListN<3>& gts_triangle_list)	// out, vertex triples, anti-clockwise
    {
        // open the file
        std::ifstream gts_file;
        gts_file.open(gts_filename.c_str());
        assert(gts_file.good());
		constexpr auto maxnum = std::numeric_limits<std::streamsize>::max();

        // read in gts-format data
        unsigned int num_verts, num_edges, num_faces;
        while (gts_file.peek() == '#' || gts_file.peek() == '!') { // skip cmnts
			gts_file.ignore(maxnum, '\n');
			assert(!gts_file.eof());
		}
        gts_file >> num_verts >> num_edges >> num_faces;
        gts_file.ignore(maxnum, '\n');
        assert(num_edges == (3*num_faces / 2));

        // get vertices
        vertex_list.reserve(num_verts);
        double vx,vy,vz;
        double nx=0,ny=0,nz=0;
        for (unsigned int vidx=0; vidx < num_verts; ++vidx) {
            nx=0;
            ny=0;
            nz=0;
            while (gts_file.peek() == '#' || gts_file.peek() == '!') {
				gts_file.ignore(maxnum, '\n');
				assert(!gts_file.eof());
			}
            std::string line;
            getline(gts_file, line);
            std::stringstream buf(line);
            buf >> vx >> vy >> vz >> nx >> ny >> nz;
            vertex_list.push_back(Vertex(Vector(vx,vy,vz), Vector(nx,ny,nz)));
        }
        //std::cout << "Added " << vertex_list.size() << " vertices." << std::endl;

        // get edge list
        gts_edge_list.reserve(num_edges);
        std::array<unsigned int,2> edge_vert;
        for (unsigned int eidx=0; eidx < num_edges; ++eidx) {
            while (gts_file.peek() == '#' || gts_file.peek() == '!') {
				gts_file.ignore(maxnum, '\n');
				assert(!gts_file.eof());
			}
            gts_file >> edge_vert[0] >> edge_vert[1];
            gts_file.ignore(maxnum, '\n');
			for (int i = 0; i < 2; i++) edge_vert[i]--;
            gts_edge_list.push_back(edge_vert);
		}

		// get triangle list
        gts_triangle_list.reserve(num_faces);
        unsigned int e[3];
        for (unsigned int tidx=0; tidx < num_faces; ++tidx) {
            while (gts_file.peek() == '#' || gts_file.peek() == '!') {
				gts_file.ignore(maxnum, '\n');
				assert(!gts_file.eof());
			}
            gts_file >> e[0] >> e[1] >> e[2];
            gts_file.ignore(maxnum, '\n');

			// Translate three edges to three vertices in order
			// Get the common vertex of each edge pair
			constexpr unsigned int perm[] = { 1, 2, 0 };
			for (int i = 0; i < 3; i++) e[i] -= 1;
			std::vector<unsigned int> intersection;
			std::array<unsigned int,3> tri;
			for (int i = 0; i < 3; i++) {
				std::array<unsigned int,2>& vl1 = gts_edge_list[e[i]];
				std::array<unsigned int,2>& vl2 = gts_edge_list[e[perm[i]]];
				intersection.clear();
				std::set_intersection(vl1.cbegin(), vl1.cend(),
									  vl2.cbegin(), vl2.cend(),
						std::inserter(intersection, intersection.begin()));
				tri[i] = *intersection.begin();
			}
            gts_triangle_list.push_back(tri);
		}
	};

	static void write_mesh_to_file(
		const std::string& gts_filename,	//!\param file to write to
		const std::vector<Vertex>& vertex_list,	//!\param in vertices
		const IdxListN<2>& edge_list,		//!\param in vertex indices low-high
		const IdxListN<3>& triangle_list)	//!\param in edge indices anti-clock
	{
        // open the file
        std::ofstream gts_file;
        gts_file.open(gts_filename.c_str());
        assert(gts_file.good());
		constexpr auto maxnum = std::numeric_limits<std::streamsize>::max();

        // write out gts-format data
        unsigned int num_verts = vertex_list.size();
		unsigned int num_edges = edge_list.size();
		unsigned int num_faces = triangle_list.size();
        assert(num_edges == (3*num_faces / 2));
		gts_file << num_verts << " " << num_edges << " " << num_faces
				 << " GtsSurface GtsFace GtsEdge GtsVertex" << std::endl;

        // write vertices
        for (const auto& v: vertex_list) {
			const Vector& n = v.get_normal();
            gts_file << v.x << " " << v.y << " " << v.z << " "
					 << n.x << " " << n.y << " " << n.z << std::endl;
        }

        // write edge list
        for (const auto& ev: edge_list) {
            gts_file << (ev[0]+1) << " " << (ev[1]+1) << std::endl;;
		}

		// write triangle list
        for (const auto& tv: triangle_list) {
            gts_file << (tv[0]+1) << " " << (tv[1]+1) << " " << (tv[2]+1)
					 << std::endl;
		}
	};
};

#endif // _BEEP_GTS_UTILS_H_
