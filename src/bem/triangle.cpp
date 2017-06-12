/*
* triangle.cpp
*
*  Created on: 21 Jul 2010
*      Author: david
*/

#include "triangle.h"

  Vector get_field(const BasicTriangle& tri, const Vector& fvals, const Vector& hvals, double u, double v)
    {
#if 1
        //std::cout << "get_field @ (u,v) = (" << u << "," << v << ")" << "\n";
        //Vector E = f_at_verts_to_E(tri.get_v1(), tri.get_v2(), tri.get_v3(), fvals.x, fvals.y, fvals.z); 
        const Vector& v1 = tri.get_v1();
        const Vector& v2 = tri.get_v2();
        const Vector& v3 = tri.get_v3();
        double f1 = fvals.x;
        double f2 = fvals.y;
        double f3 = fvals.z;
        double w = 1.0 - u - v;

        Vector axis1 = (v2 - v1).normalised();
        Vector axis2 = (v3 - v1).normalised();
        double df_axis1 = (f2 - f1) / (v2-v1).length();
        double df_axis2 = (f3 - f1) / (v3-v1).length();
        //return - ((axis1*df_axis1) + (axis2*df_axis2));
        
        //Vector pt = tri.get_point(u,v);
        //double f_interp = w*fvals[0] + u*fvals[1] + v*fvals[2];
        //double f_delta_u = ((w-eps/2.0)*fvals[0] + (u+eps)*fvals[1] + (v-eps/2.0)*fvals[2]) - f_interp;
        //double f_delta_v = ((w-eps/2.0)*fvals[0] + (u-eps/2.0)*fvals[1] + (v+eps)*fvals[2]) - f_interp;
        //double f_delta_u = (fvals[1] - fvals[2]*0.5 - fvals[0]*0.5) / (tri.get_v2() - tri.get_v1()*0.5 - tri.get_v3()*0.5).length();
        //double f_delta_v = (fvals[2] - fvals[0]*0.5 - fvals[1]*0.5) / (tri.get_v3() - tri.get_v1()*0.5 - tri.get_v2()*0.5).length();
        
        //Vector dU  = tri.get_dU(u,v,w);
        //dU.normalise();
        //std::cout << "dU (numer/analytic): " << dU << " " << tri.get_dU(u,v,w) << "\n";
        //Vector dV = (tri.get_point(u-eps/2.0,v+eps) - pt);
        //Vector dV = tri.get_dV(u,v,w);
        //dV.normalise();
        //std::cout << "dV (numer/analytic): " << dV << " " << tri.get_dV(u,v,w) << "\n";
        //Vector n = dU.cross(dV).normalised();
        
        //if (n.dot(tri.get_planar_normal()) < 0)
        //{
        //    n = -n;
            //std::cout << tri.get_v1() << " " << tri.get_v2() << " " << tri.get_v3() << "\n";
        //}

        //Vector hn = w*hvals[0]*tri.get_n1() + u*hvals[1]*tri.get_n2() + v*hvals[2]*tri.get_n3();
        //double h = hn.length();
        //hn.normalise();
        double h = w*hvals[0] + u*hvals[1] + v*hvals[2];
        Vector hn = w*tri.get_n1() + u*tri.get_n2() + v*tri.get_n3();
        hn.normalise();
        
        //Vector n_interp = w*tri.get_n1() + u*tri.get_n2() + v*tri.get_n3();
        //n_interp.normalise();

        double *a_data = new double[9];
        double *f_data = new double[3];
        gsl_matrix_view m = gsl_matrix_view_array(a_data, 3, 3);
        gsl_vector_view f = gsl_vector_view_array(f_data, 3);

        gsl_matrix_set(&m.matrix, 0, 0, axis1.x);
        gsl_matrix_set(&m.matrix, 0, 1, axis1.y);
        gsl_matrix_set(&m.matrix, 0, 2, axis1.z);
        gsl_matrix_set(&m.matrix, 1, 0, axis2.x);
        gsl_matrix_set(&m.matrix, 1, 1, axis2.y);
        gsl_matrix_set(&m.matrix, 1, 2, axis2.z);
        gsl_matrix_set(&m.matrix, 2, 0, hn.x);
        gsl_matrix_set(&m.matrix, 2, 1, hn.y);
        gsl_matrix_set(&m.matrix, 2, 2, hn.z);

        gsl_vector_set(&f.vector, 0, df_axis1);
        gsl_vector_set(&f.vector, 1, df_axis2);
        gsl_vector_set(&f.vector, 2, h);

        // do the SVD
        gsl_vector *solution = gsl_vector_alloc(3);
        gsl_vector *work = gsl_vector_alloc(3);
        gsl_matrix *V = gsl_matrix_alloc(3,3);
        gsl_vector *S = gsl_vector_alloc(3);

        int retval = gsl_linalg_SV_decomp(&(m.matrix), V, S, work);
        if (retval != 0) { std::cerr << "SVD_decomp retval: " << retval << std::endl; }
        assert(retval == 0);
        retval = gsl_linalg_SV_solve(&m.matrix, V, S, &f.vector, solution);
        if (retval != 0) { std::cerr << "SVD_solve retval: " << retval << std::endl; }
        assert(retval == 0);
    
        Vector E;
        E.x = gsl_vector_get(solution, 0);
        E.y = gsl_vector_get(solution, 1);
        E.z = gsl_vector_get(solution, 2);

        //std::cout << "E: " << E << " (h.n= " << h_interp*n << ")" << "\n";

        // delete allocated memory
        gsl_vector_free(solution);
        gsl_vector_free(work);
        gsl_matrix_free(V);
        gsl_vector_free(S);
        delete[] a_data;
        delete[] f_data;
        
        return E;
#else
        Vector E = f_at_verts_to_E(tri.get_v1(), tri.get_v2(), tri.get_v3(), fvals.x, fvals.y, fvals.z); 
        
        double w = 1.0 - u - v;
        Vector n = tri.get_normal(u,v);
        double h_interp = w*hvals[0] + u*hvals[1] + v*hvals[2];
        E += h_interp*n;
        
        return E;
#endif
    }

     Vector ionic_from_E(const Vector& Eo,
                       const Vector& Ei,
                       const Vector& norm,
                       double kappa,
                       double f,
                       double epsilon_int,
                       double epsilon_ext)
    {
        
        double o_i = epsilon_ext - epsilon_int;
        Vector ionic_pressure = -0.5*epsilon_ext*f*f*kappa*kappa*norm;
        return ionic_pressure;
    }

    Vector MST_from_E(const Vector& E,
                       const Vector& norm,
                       double dielectric)
    {
        const double En = E.dot(norm);
        double E2 = E.dot(E);
        Vector mst(0,0,0);
        mst.x = dielectric*(E.x*En - 0.5*E2*norm.x);
        mst.y = dielectric*(E.y*En - 0.5*E2*norm.y);
        mst.z = dielectric*(E.z*En - 0.5*E2*norm.z);
        return mst;
    }

    // these two functions should give identical values...
    Vector dbf_from_MST(const Vector& Eo,
                        const Vector& Ei,
                        const Vector& norm,
                        double epsilon_int,
                        double epsilon_ext)
    {
        return MST_from_E(Eo, norm, epsilon_ext) - MST_from_E(Ei, norm, epsilon_int);
    }
     Vector dbf_from_E(const Vector& Eo,
                       const Vector& Ei,
                       const Vector& norm,
                       double epsilon_int,
                       double epsilon_ext)
    {
        double o_i = epsilon_ext - epsilon_int;
        Vector dbf = -0.5*o_i*(Eo.dot(Ei))*norm;
        return dbf;
    }

    Vector dbf_from_delta(const Vector& Eo,
                          const Vector& Ei,
                          const Vector& norm,
                          double epsilon_int,
                          double epsilon_ext)
    {
        double prefactor = -0.5*(epsilon_ext - epsilon_int);
        double zo = Eo.dot(norm);
        double zi = Ei.dot(norm);
        double brackets = Ei.dot(Ei) + zi*(zo-zi) + (zo-zi)*(zo-zi)/3.0;
        return norm * prefactor * brackets;
    }

    void calculate_force_components(const BasicTriangle& tri,
                                     Vector fvals,
                                     Vector hvals,
                                     double epsilon_int,
                                     double epsilon_ext,
                                     double kappa,
                                     const QuadratureRule& rule,
                                     unsigned int subdivides,
                                     KahanVector& MST_external,
                                     KahanVector& MST_internal,
                                     KahanVector& dbf,
                                     KahanVector& ionic)
    {
        // Generate a set of parametric quadrature points at which to evaluate force
        // component over the triangle.  Use double precision.
        typedef std::vector< QuadPointT<double> > DblQuadList;
        
        DblQuadList parametric_quadrature_points;
        BasicTriangle parametric_triangle(Vector(1,0,0),Vector(0,1,0),Vector(0,0,1));
        parametric_triangle.get_quad_points(parametric_quadrature_points, rule, subdivides);
        
        double dielectric_ratio = epsilon_ext / epsilon_int;

        for (DblQuadList::iterator it=parametric_quadrature_points.begin(), end=parametric_quadrature_points.end(); it != end; ++it)
        {
            // parametric coordinates on the triangle
            const Vector& pt = it->pt();
            const double& u = pt.x;
            const double& v = pt.y;
            const double& w = pt.z;
            //double weight = it->weight() * area_weighting;
            
            Vector norm = tri.get_normal(u, v);
            double weight = tri.get_area() * it->weight() / parametric_triangle.get_area();
            
            Vector Eout = get_field(tri, fvals, hvals, u, v);
            Vector Ein = Eout + (Eout.dot(norm)*(dielectric_ratio - 1.0))*norm;
            
            //Vector Ein = get_field(tri, fvals, hvals*dielectric_ratio, u, v);
            //Vector Eout = Ein - (Ein.dot(norm)*(1.0 - 1.0/dielectric_ratio))*norm;
            
            MST_external += (MST_from_E(Eout, norm, epsilon_ext) * weight);
            MST_internal += (MST_from_E(Ein, norm, epsilon_int) * weight);
            //double delta = solve_delta(Eout, Ein, norm);
            //dbf += dbf_from_delta(Eout, Ein, norm, epsilon_int, epsilon_ext) * weight;
            dbf += dbf_from_E(Eout, Ein, norm, epsilon_int, epsilon_ext) * weight;
            ionic += ionic_from_E(Eout, Ein, norm, kappa, fvals.x*w + fvals.y*u + fvals.z*v, epsilon_int, epsilon_ext) * weight;
        }
        
        return;
    }

     Vector f_at_verts_to_E(const Vector& vertex1, const Vector& vertex2, const Vector& vertex3,
                            double vf1, double vf2, double vf3)
    {

        const Vector& v1 = vertex1;
        const Vector& v2 = vertex2;
        const Vector& v3 = vertex3;
        const double& f1 = vf1;
        const double& f2 = vf2;
        const double& f3 = vf3;

        Vector axis1 = (v2 - v1).normalised();
        //assert(fabs(normal.dot(axis1)) < 1e-6);
        Vector normal = axis1.cross(v3-v1).normalised();
        Vector axis2 = normal.cross(axis1).normalised();
        double df_axis1 = (f2-f1) / (v2-v1).length();
        double df_axis2 = ((f3-f1) - axis1.dot(v3-v1)*df_axis1) / axis2.dot(v3-v1);
        return - ((axis1*df_axis1) + (axis2*df_axis2));

//         Vector E(0,0,0);
//         // permutations of vector triplets
//         const Vector* vertices[3] = {&vertex1, &vertex2, &vertex3};
//         const double* values[3]   = {&vf1, &vf2, &vf3};
//         
//         gsl_permutation *p = gsl_permutation_alloc(3);
//         gsl_permutation_init(p);
//         int ctr=0;
//         do
//         {
//             int idx1 = gsl_permutation_get(p,0);
//             int idx2 = gsl_permutation_get(p,1);
//             int idx3 = gsl_permutation_get(p,2);
// 
//             const Vector& v1 = *(vertices[idx1]);
//             const Vector& v2 = *(vertices[idx2]);
//             const Vector& v3 = *(vertices[idx3]);
//             const double& f1 = *(values[idx1]);
//             const double& f2 = *(values[idx2]);
//             const double& f3 = *(values[idx3]);
// 
//             Vector axis1 = (v2 - v1).normalised();
//             //assert(fabs(normal.dot(axis1)) < 1e-6);
//             Vector normal = axis1.cross(v3-v1).normalised();
//             Vector axis2 = normal.cross(axis1).normalised();
//             double df_axis1 = (f2-f1) / (v2-v1).length();
//             double df_axis2 = ((f3-f1) - axis1.dot(v3-v1)*df_axis1) / axis2.dot(v3-v1);
//             E = E - (axis1*df_axis1) - (axis2*df_axis2);
// 
//             ctr++;
//         }
//         while (gsl_permutation_next(p) == GSL_SUCCESS);
//         gsl_permutation_free (p);

//         // expecting 6 permutations
//         assert(ctr == 6);
// 
//         E = E / ctr;
// 
//     /*    std::cout << (v2-v1).dot(E) << " " << (f2-f1) << std::endl;
//         std::cout << (v3-v1).dot(E) << " " << (f3-f1) << std::endl;
//         std::cout << "h,E,f1,f2,f3,err_ax1,err_ax2: " << h << " " << E << " " << f1 << " " << f2 << " " << f3 << " " << fabs(f2 - f1 + (v2-v1).dot(E)) << " " << fabs((f3 - f1) - -(v3-v1).dot(E)) << std::endl;
//         assert(fabs((f2 - f1) - -(v2-v1).dot(E)) < 1e-6);
//         assert(fabs((f3 - f1) - -(v3-v1).dot(E)) < 1e-6);*/
//         //assert(E.dot(normal) >= 0);
// 
//         return E;

    }
