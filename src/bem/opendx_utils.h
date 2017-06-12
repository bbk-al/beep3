#ifndef __OPENDX_UTILS_H_
#define __OPENDX_UTILS_H_

class GridParms
{
    
public:
        
    // constructor for GridParms
    GridParms()
    {
        this->x_pts = 257;
        this->y_pts = 257;
        this->z_pts = 257;
        
        this->x_angstroms = 0.4f*(x_pts - 1);
        this->y_angstroms = 0.4f*(y_pts - 1);
        this->z_angstroms = 0.4f*(z_pts - 1);
        
        this->origin_x = -x_angstroms / 2.0;
        this->origin_y = -y_angstroms / 2.0;
        this->origin_z = -z_angstroms / 2.0;
        
        opendx_suffix = "attribute \"dep\" string \"positions\"\n"
                "object \"regular positions regular connections\" class field\n"
                "component \"positions\" value 1\n"
                "component \"connections\" value 2\n"
                "component \"data\" value 3\n";
        
    };    
    
    Vector
    GridIdxToPos(unsigned int x, unsigned int y, unsigned int z) const
    {
        
        const double xpos = (x*this->deltax()) + origin_x;
        const double ypos = (y*this->deltay()) + origin_y;
        const double zpos = (z*this->deltaz()) + origin_z;

        return Vector(xpos,ypos,zpos);

    };

    std::string
    OpenDX_Preamble() const
    {
        std::ostringstream buf;
        buf << std::scientific; // get scientific notation
        buf <<  "# Data from BEM\n"
                "#\n"
                "# POTENTIAL (kT/e)\n"
                "#\n"
                "object 1 class gridpositions counts " <<  this->x_pts << " " << this->y_pts << " " << this->z_pts << "\n"
                "origin " << this->origin_x << " " << this->origin_y << " " << this->origin_z << "\n"
                "delta " << this->deltax() << " 0.000000e+00 0.000000e+00\n"
                "delta 0.000000e+00 " << this->deltay() << " 0.000000e+00\n"
                "delta 0.000000e+00 0.000000e+00 " << this->deltaz() << "\n"
                "object 2 class gridconnections counts " << this->x_pts << " " << this->y_pts << " " <<  this->z_pts << "\n"
                "object 3 class array type double rank 0 items " << this->num_grid_points() << " data follows\n";
        return buf.str();
    };
    
    std::string
    OpenDX_Suffix() const
    {
        return this->opendx_suffix;
    }
    
    double deltax() const { return this->x_angstroms / double(this->x_pts - 1); };
    double deltay() const { return this->y_angstroms / double(this->y_pts - 1); };
    double deltaz() const { return this->z_angstroms / double(this->z_pts - 1); };

    unsigned int num_grid_points() const { return this->x_pts * this->y_pts * this->z_pts; };
    
    unsigned int x_pts;
    unsigned int y_pts;
    unsigned int z_pts;
    double x_angstroms;
    double y_angstroms;
    double z_angstroms;
    double origin_x;
    double origin_y;
    double origin_z;
    
    std::string opendx_suffix;
};

#endif
