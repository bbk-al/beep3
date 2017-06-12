/*
* vertex.cpp
*
*  Created on: 21 Jul 2010
*      Author: david
*/

#include "vertex.h"
#include "triangle.h"
#include <sstream>
#include <string>
#include <gsl/gsl_linalg.h>
#include <vector>
#include "../common/math_vector.h"

static Vector centroid_from_points(const std::vector<Vector>& points)
{
    Vector centroid(0,0,0);
    for (std::vector<Vector>::const_iterator it=points.begin(),end=points.end();
        it != end;
        ++it)
    {
        centroid += *it;
    }
    centroid = centroid / points.size();
    return centroid;
}

static Vector orthogonal_distance_regression(const std::vector<Vector>& points, const Vector& normal_approx)
{

    Vector centroid = centroid_from_points(points);

    // Orthogonal Distance Regression for points in BasicNodePatch
    // This will give us minimal plane for points -- use this to get a smoother
    // set of patch normals from which to estimate the solid angle subtended by
    // the node point.

    int dim = points.size();
    gsl_matrix *m = gsl_matrix_alloc(dim, 3);
    //gsl_vector_view vector_b = gsl_vector_view_array(b_data, dim);
    //gsl_vector *solution = gsl_vector_alloc(2);
    gsl_vector *work = gsl_vector_alloc(3);
    gsl_matrix *V = gsl_matrix_alloc(3, 3);
    gsl_vector *S = gsl_vector_alloc(3);

    // set values of m
    unsigned int ctr=0;
    for (std::vector<Vector>::const_iterator it=points.begin(),end=points.end();
        it != end;
        ++it)
    {
        const Vector& v = *it;
        gsl_matrix_set(m, ctr, 0, v.x - centroid.x);
        gsl_matrix_set(m, ctr, 1, v.y - centroid.y);
        gsl_matrix_set(m, ctr, 2, v.z - centroid.z);
        ctr++;
    }
    assert(ctr == points.size());

    // get singular value decomposition of M
    int retval = gsl_linalg_SV_decomp(m,V,S,work);
    assert(retval == 0);

    // normal for plane
    Vector nfit(0,0,0);
    nfit.x = gsl_matrix_get(V, 0, 2);
    nfit.y = gsl_matrix_get(V, 1, 2);
    nfit.z = gsl_matrix_get(V, 2, 2);
    nfit.normalise();

    if (nfit.dot(normal_approx) < 0.0) { nfit = nfit * -1; }

    // extract eigenvector corresponding to smallest eigenvalue
    //retval = gsl_linalg_SV_solve(&(matrix_a.matrix),V,S,&(vector_b.vector),solution);
    //assert(retval == 0);

    gsl_matrix_free(m);
    gsl_vector_free(work);
    gsl_matrix_free(V);
    gsl_vector_free(S);

    return nfit;
}

void Vertex::set_normal(const std::vector<Triangle>& all_triangles)
{
   
#if 0
    // init rough guess
    normal = Vector(0.0,0.0,0.0);
    std::vector<Vector> pts;
    for(size_t i=0; i < triangle_indices.size(); ++i )
    {
        size_t t_idx = triangle_indices[i];
        const Triangle& tri = all_triangles[t_idx];
        normal = normal + tri.get_planar_normal()*tri.get_planar_area();
        pts.push_back(tri.get_v1());
        pts.push_back(tri.get_v2());
        pts.push_back(tri.get_v3());
    }
    normal.normalise();

    // remove central node vertex
    std::vector<Vector> pts2;
    for (std::vector<Vector>::const_iterator it=pts.begin(), end=pts.end(); it != end; ++it)
    {
        if (vectors_are_equal(*it, *this) == false)
        {
            pts2.push_back(*it);
        }
    }

    // now guess better using orthogonal distance regression of neighbouring points
    normal = orthogonal_distance_regression(pts2, normal);
#else

    // init rough guess
    normal = Vector(0.0, 0.0, 0.0);
    for(size_t i=0; i < triangle_indices.size(); ++i )
    {
        size_t t_idx = triangle_indices[i];
        const Triangle& tri = all_triangles[t_idx];
        
        // ignore the influence of very small triangles
        if (tri.get_planar_area() < 1e-12) { continue; } 
        
        normal = normal + tri.get_planar_normal()*tri.get_planar_area();
    }
    
    if (normal.length2() < 1e-12)
    {
        // oh dear this vertex looks highly suspicious
        std::cerr << "Bad Vertex: " << *this << std::endl;
        return;
    }
    
    normal.normalise();

#endif

    return;
}

void Vertex::reorder(const std::vector<Triangle>& all_triangles)
{

    // reorders the triangular elements of the NodePatch
    for(size_t i=0; i < triangle_indices.size()-1; ++i )
    {
        size_t t_idx = triangle_indices[i];
        const Triangle& tri = all_triangles[t_idx];

        Vector a;
        Vector b;
        tri.get_anti_clockwise_vertices_from(*this, a, b);

        for (size_t j=i+1; j < triangle_indices.size(); ++j)
        {

            size_t t2_idx = triangle_indices[j];
            const Triangle& tri_2 = all_triangles[t2_idx];
            Vector c;
            Vector d;
            tri_2.get_anti_clockwise_vertices_from(*this, c, d);

            if (vectors_are_equal(b,c))
            {
                std::swap(triangle_indices[i+1], triangle_indices[j]);
                break;
            }

        }
    }

    return;

}

std::string Vertex::kinemage_edges(const std::vector<Triangle>& all_triangles) const
{
    const Vertex& v = *this;
    
    std::ostringstream buf;
    buf << "{}P ";
    
    for(size_t i=0; i < triangle_indices.size(); ++i )
    {
        
        size_t t_idx = triangle_indices[i];
        const Triangle& tri = all_triangles[t_idx];
    
        Vector a;
        Vector b;
        try {
            
            tri.get_anti_clockwise_vertices_from(v, a, b);

            Vector mid_a = (v + a) / 2.0;
            Vector mid_b = (v + b) / 2.0;
            Vector centre = tri.get_planar_centroid();

            if (i==0) {
                buf << mid_a.x << " " << mid_a.y << " " << mid_a.z << " ";
            }
            buf << centre.x << " " << centre.y << " " << centre.z << " ";
            buf << mid_b.x << " " << mid_b.y << " " << mid_b.z << " ";
        }
        catch (Vertex::BadVertex &bv) {
            std::cout << "Vertex " << v << " not a member of Triangle " << tri << ".  Puzzling ... " << std::endl;
        }
        
    }
    buf << "\n";
    return buf.str();
}


