/*
* generalized_born.h
*
*  Created on: 21 Jul 2010
*      Author: david
*/
#ifndef _GENERALIZED_BORN_H_
#define _GENERALIZED_BORN_H_

#include <math.h>
#include <cassert>
#include <string>
#include <vector>
#include "../bem/triangle.h"
#include "../bem/vertex.h"
#include "../bem/gts_utils.h"
#include "../common/charge.h"
#include "../bem/quad_point.h"
#include "../bem/png1.h"
#include "../bem/constants.h"

class GeneralizedBorn
{

public:

    GeneralizedBorn(const std::string gts_filename, const std::string& charges_filename)
    {
        // only need these whilst reading in mesh
        std::vector<Vertex> vertices;
        std::vector<Triangle> triangles;

        Charge::read_charges_from_file(charges_filename, charges, true); // do not skip zero charges- need to know number of atoms
        GTS_Surface::read_mesh_from_file(gts_filename, vertices, triangles);

        // Now loop over the triangles and convert them into Bezier triangles,
        // then extract quadrature points from them.
        for (std::vector<Triangle>::const_iterator tri_it=triangles.begin(), tri_end=triangles.end();
             tri_it != tri_end;
             ++tri_it)
        {
            // create bezier triangle
            try {
                PNG1_Triangle png(*tri_it);
                png.get_quad_points(surface_quads, BasicTriangle::gauss_legendre_pts_4());
            }
            catch (BadQuadPoint)
            {
                std::cerr << "WARNING: Bad QuadPoint generated on Triangle: " << *tri_it << " (skipping triangle)" << std::endl;
            }
        }

        born_radii = new double[charges.size()];
        assert(born_radii != NULL);
        
        double total_area = 0.0;
        std::cout << "Got " << triangles.size() << " faces (" << vertices.size() << " vertices) and " << charges.size() << " charges.  Generated " << surface_quads.size() << " quadrature points." << std::endl;
        //std::cout << "@kinemage\n@vectorlist {quadrature_points}\n";
        for (std::vector<QuadPoint>::iterator it=surface_quads.begin(), end=surface_quads.end(); it != end; ++it)
        {
            const QuadPoint& qp = *it;
            //Vector a = qp.pt();
            //Vector b = a + qp.normal() * qp.weight();
            total_area += qp.weight();
            //std::cout << "{}P " << a.x << " " << a.y << " " << a.z << " " << b.x << " " << b.y << " " << b.z << "\n";
        }
        
        std::cout << "Total surface area: " << total_area << std::endl;
    }

    ~GeneralizedBorn() {
        delete[] born_radii;
    }

    void calculate_born_radii()
    {
        unsigned int charge_ctr=0;
        for (std::vector<Charge>::const_iterator ch_it=charges.begin(), ch_end=charges.end(); ch_it != ch_end; ++ch_it)
        {
            double surf_int = 0;
            const Vector& x = ch_it->position();

            if (ch_it->get_charge() != 0)
            {
                for (std::vector<QuadPoint>::const_iterator qp_it=surface_quads.begin(), qp_end=surface_quads.end();
                    qp_it != qp_end;
                    ++qp_it)
                {
                    const QuadPoint& qp = *qp_it;
                    Vector r(qp.pt());
                    r -= x;
                    double r2 = r.length2();
                    double r4 = r2*r2;
                    surf_int += qp.weight() * r.dot(qp.normal()) / r4;
                }
                born_radii[charge_ctr++] = 1.0 / (surf_int * ONE_OVER_4PI);
            }
            else
            {
                born_radii[charge_ctr++] = 0;
            }
            
        }

    }
    
    void calculate_born_radii_GBR6()
    {
        unsigned int charge_ctr=0;
        for (std::vector<Charge>::const_iterator ch_it=charges.begin(), ch_end=charges.end(); ch_it != ch_end; ++ch_it)
        {
            double surf_int = 0;
            const Vector& x = ch_it->position();

            if (ch_it->get_charge() != 0)
            {

                for (std::vector<QuadPoint>::const_iterator qp_it=surface_quads.begin(), qp_end=surface_quads.end();
                    qp_it != qp_end;
                    ++qp_it)
                {
                    const QuadPoint& qp = *qp_it;
                    Vector r(qp.pt());
                    r -= x;
                    double r2 = r.length2();
                    double r6 = r2*r2*r2;
                    surf_int += qp.weight() * r.dot(qp.normal()) / r6;
                    
                }
                born_radii[charge_ctr++] = pow(3.0 * surf_int * ONE_OVER_4PI, -1.0/3.0);
            }
            else
            {
                born_radii[charge_ctr++] = 0;
            }
            
        }

    }
    
    double lpb_gbr6_scaling(double epsilon_int, double epsilon_ext) const
    {
        unsigned int num_atoms = charges.size();
        double net_charge = 0;
        for (std::vector<Charge>::const_iterator ch_it=charges.begin(), ch_end=charges.end(); ch_it != ch_end; ++ch_it)
        {
            net_charge += ch_it->get_charge();
        }

        double ratio = epsilon_int / epsilon_ext;
        double charge_term = pow(fabs(net_charge), 0.65);
        double A = (-1.63e-3*charge_term) + (2.18e-6*num_atoms) + 1.016;
        double B = ( 3.31e-2*charge_term) - (4.77e-5*num_atoms) + 0.683;
        double numerator = A + (2.0*B*ratio);
        double denominator = 1.0 + (2.0*ratio);
        return numerator / denominator;
    }

    double evaluate_energy(double epsilon_int, double epsilon_ext, double F=4.) const
    {
        double GB_energy = 0;
        unsigned int ii=0, jj=0;
        for (unsigned int ii=0; ii < charges.size(); ++ii)
        {
            const Charge& chi = charges[ii];
            if (chi.get_charge() == 0) continue;
            double Ri = born_radii[ii];

            for (unsigned int jj=ii; jj < charges.size(); ++jj)
            {
                const Charge& chj = charges[jj];
                if (chj.get_charge() == 0) continue;
                double Rj = born_radii[jj];
                //std::cout << "ri,rj: " << Ri << " " << Rj << std::endl;

                double numerator = chi.get_charge() * chj.get_charge();
                double rij2 = (chj - chi).length2();
                double RiRj = Ri * Rj;
                double exp_term = exp(-rij2/(F*RiRj));
                double denom = sqrt(rij2 + RiRj*exp_term);

                double multiplic = (ii == jj) ? 1.0 : 2.0;

                GB_energy += multiplic * numerator / denom;

            }
        }

        GB_energy *= -0.5*ONE_OVER_4PI*(1./epsilon_int - 1./epsilon_ext);

        return GB_energy;

    }

private:

    std::vector<Charge> charges;
    std::vector<QuadPoint> surface_quads;
    double* born_radii;

};

#endif  // _GENERALIZED_BORN_H_
