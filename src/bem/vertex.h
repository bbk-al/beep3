/*
* vertex.h
*
*  Created on: 21 Jul 2010
*      Author: david
*/

#ifndef VERTEX_H_
#define VERTEX_H_

#include "../common/math_vector.h"
#include <vector>
#include <string>
#include <exception>

// fwd decls.
class Triangle;
class Mesh;

class Vertex : public Vector
{
    
public:

    friend class Mesh;
    friend class Triangle;

    Vertex() : Vector(0,0,0), normal(0,0,0) {}
    Vertex(const Vector& v) : Vector(v), normal(0,0,0) {}
    Vertex(const Vector& v, const Vector& vn) : Vector(v), normal(vn) {}
    
    // copy ctor
    Vertex(const Vertex& other) : Vector(static_cast<const Vector&>(other)), normal(other.normal) {
        triangle_indices.insert(triangle_indices.begin(), other.triangle_indices.begin(), other.triangle_indices.end());
    }
    
    virtual ~Vertex() {}

    void add_triangle(unsigned int tri_idx)
    {
        triangle_indices.push_back(tri_idx);
    }

    void set_normal(const std::vector<Triangle>&);

    inline const Vector& get_vertex() const { return *this; }
    inline const Vector& get_normal() const { return normal; }
    inline const std::vector<unsigned int>& get_triangle_indices() const { return triangle_indices; }

    const Vector& operator() () const { return *this; }

    // Exception definition for when user passes a bad vertex
    // intended for use with get_anti_clockwise_vertices_from()
    class BadVertex : public std::exception
    {
    public:
        BadVertex() : std::exception() {}

    };

    void reorder(const std::vector<Triangle>&);
    std::string kinemage_edges(const std::vector<Triangle>& all_triangles) const;
    
protected:

    std::vector<unsigned int> triangle_indices;
    Vector normal;
};

#endif /* VERTEX_H_ */
