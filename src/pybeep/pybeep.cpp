
#include <boost/python.hpp>
#include "../common/math_vector.h"
#include "../bem/mesh.h"
#include "../bem/node_patch.h"
#include "../bem/opendx_utils.h"
#include "../fmm/octree.h"
#include "../beep/beep.h"
#include <sstream>
#include "../fmm/interaction_list.h"

using namespace boost::python;
template<class T>
class BasicOctree : public Octree<Node<T>, T>
{
    typedef Octree<Node<T>, T> super;

public:

    BasicOctree(unsigned int max_items, const Vector& centre, double edge_length) : super(max_items, centre, edge_length) {
        num_octree_nodes=0;
        num_explicit_calcs=0;
        max_octree_depth=0;
        num_fmm_evals=0;
        num_interaction_list_moves=0;
    }

    // synonym for add -- boost::python gets confused by overloading
    inline void insert(T& thing) { super::add(&thing); }
    inline void finalize() {
        super::build_neighbourhoods();
        super::optimize_neighbourhoods();
        super::remove_empty_nodes();
        
        num_octree_nodes = super::get_total_num_nodes();
        num_explicit_calcs = super::calc_neighbourhood_interacts();
        max_octree_depth = super::get_bottom_level() - super::get_top_level();
        num_fmm_evals = (max_octree_depth > 1) ? super::get_total_items() : 0;
        num_interaction_list_moves = calc_total_interaction_list_moves();
      
    }

    int check_for_duplicate(const Vector& xyz) const
    {
        const typename super::NodeT::ContentList& contents = super::get_everything();
        for (int ii=0; ii < super::size(); ++ii)
        {
            Vector pt = static_cast<Vector&>(*(contents[ii]));
            if ((pt - xyz).length2() < 1e-9)
            {
                return ii;
            }
        }
        return -1;
    }
    
    size_t calc_total_interaction_list_moves() const
    {
        // loop over all nodes, don't bother with leaf nodes or sub-leaf nodes
        // add up the interactions between children of these nodes
        size_t ilist=0;
        
        // in the FMM near-field.
        for (int level=super::get_top_level(); level <= super::get_bottom_level() && level < super::max_depth; ++level)
        {
            for (typename super::NodeList::const_iterator node_it=super::get_node_list(level).begin(), node_end=super::get_node_list(level).end(); 
                 node_it != node_end; ++node_it)
            {
                const typename super::NodeT& node = *(node_it->second);
                if (node.isLeaf() == false && node.isShadow() == false) 
                {
                    // get interaction list for this neighbour_id
                    for (unsigned short ctr=0; ctr < ILIST_SIZE; ++ctr)
                    {
                        const INTERACTION& inter = interactions[ctr];
                        if (inter.neighbour == 13) { continue; }
                        try {
                     
                            OctreeIndexer idxer1 = node.get_idx().get_neighbour_idxer(inter.neighbour).get_child_idx(inter.targ);
                            const typename super::NodeT& ch1 = super::get_node(idxer1);
                            OctreeIndexer idxer2 = node.get_idx().get_child_idx(inter.src);
                            const typename super::NodeT& ch2 = super::get_node(idxer2);
                            ++ilist;
                        }
                        catch (...) { }
                    }
                }
            }
        }

        return ilist;   
    }
        
    void rejig(unsigned int new_nb_size)
    {
        long start_time = myclock();
        super::rejig(new_nb_size);
        //std::cerr << "Tree re-jig took: " << (myclock() - start_time ) / 1000. << " ms" << std::endl;
        
    }
    
    std::string spew() const
    {
        std::ostringstream buf;
        for (unsigned int level=super::get_top_level(); level <= super::get_bottom_level(); ++level)
        {
            const typename super::NodeList& node_list = super::get_node_list(level);
            for (typename super::NodeList::const_iterator it=node_list.begin(), end=node_list.end(); it != end; ++it)
            {

                const OctreeIndexer& idxer = it->first;
                buf << level << " " << idxer.get_x_idx() << " " << idxer.get_y_idx() << " " << idxer.get_z_idx() << std::endl;
            }
        }
        return buf.str();
    }

    size_t num_octree_nodes;
    size_t num_explicit_calcs;
    size_t max_octree_depth;
    size_t num_fmm_evals;
    size_t num_interaction_list_moves;
    
};

BOOST_PYTHON_MODULE(libBEEP)
{
    class_<Quaternion>("Quaternion", init<double, double, double, double>())
        .def(init<const Quaternion&>())
        .def_readwrite("a", &Quaternion::a)
        .def_readwrite("b", &Quaternion::b)
        .def_readwrite("c", &Quaternion::c)
        .def_readwrite("d", &Quaternion::d)
        .def_readwrite("x", &Quaternion::a) // alternative notation for a-d
        .def_readwrite("y", &Quaternion::b)
        .def_readwrite("z", &Quaternion::c)
        .def_readwrite("w", &Quaternion::d)
        .def("__str__", &Quaternion::str)
        .def("inverse", &Quaternion::inverse);

    class_<Vector>("_Vector", init<double, double, double>())
        .def(init<const Vector&>())
        .def(init<const __Vector<float>&>())
        .def(init<const __Vector<double>&>())
        .def_readwrite("x", &Vector::x)
        .def_readwrite("y", &Vector::y)
        .def_readwrite("z", &Vector::z)
        .def("length",&Vector::length)
        .def("length2",&Vector::length2)
        .def("dot",&Vector::dot)
        .def("cross",&Vector::cross)
        .def("__normalised",&Vector::normalised)
        .def("normalise", &Vector::normalise)
        .def("__add__",&Vector::__py_add)
        .def("__sub__",&Vector::__py_subtract)
        .def("__mul__",&Vector::operator*)
        .def("__div__",&Vector::operator/)
        .def("__str__", &Vector::__str__)
        .def("apply_rotation", &Vector::apply_quaternion)
        .def("change_coordinate_frame", &Vector::py_change_coordinate_frame);

    class_<Charge>("_Charge", init<const Vector&, double, double>())
        .def(init<const Vector&, double, double>())
        .def_readwrite("charge", &Charge::charge)
        .def_readwrite("radius", &Charge::radius)
        .def_readwrite("x", &Charge::x)
        .def_readwrite("y", &Charge::y)
        .def_readwrite("z", &Charge::z)
        .def("position", &Charge::py_get_position, return_value_policy<copy_const_reference>())
        .def("_change_coordinate_frame", &Charge::py_change_coordinate_frame);
    
    class_<QuadPoint>("QuadPoint", init<const QuadPoint&>())
        .def("pt", &QuadPoint::__py_pt, return_value_policy<copy_const_reference>())
        .def("normal", &QuadPoint::__py_normal, return_value_policy<copy_const_reference>())
        .def("weight", &QuadPoint::__py_weight, return_value_policy<copy_const_reference>())
        .def("__str__", &QuadPoint::str);

    class_<BasicNodePatch>("NodePatch", init<const BasicNodePatch&>())
        .def("vector", &BasicNodePatch::py_vector, return_value_policy<copy_non_const_reference>())
        .def("node", &BasicNodePatch::get_node, return_value_policy<copy_const_reference>())
        .def("centroid", &BasicNodePatch::get_centroid, return_value_policy<copy_const_reference>())
        .def("normal", &BasicNodePatch::get_normal, return_value_policy<copy_const_reference>())
        .def("alt_normal", &BasicNodePatch::get_alt_normal, return_value_policy<copy_const_reference>())
        .def("planar_area", &BasicNodePatch::get_planar_area)
        .def("bezier_area", &BasicNodePatch::get_bezier_area)
        //.def("num_quad_points", &BasicNodePatch::get_num_quad_points)
        //.def("get_quad_point", &BasicNodePatch::get_quad_point, return_value_policy<copy_const_reference>())
        //.def("num_qual_points", &BasicNodePatch::get_num_qual_points)
        //.def("get_qual_point", &BasicNodePatch::get_qual_point, return_value_policy<copy_const_reference>())
        .def_readwrite("gc", &BasicNodePatch::gc)
        .def_readwrite("f", &BasicNodePatch::f)
        .def_readwrite("h", &BasicNodePatch::h)
        .def_readwrite("energy_coefficient_f", &BasicNodePatch::energy_coefficient_f)
        .def_readwrite("energy_coefficient_h", &BasicNodePatch::energy_coefficient_h)
        .def_readonly("force_coefficient_f", &BasicNodePatch::force_coefficient_f)
        .def_readonly("force_coefficient_h", &BasicNodePatch::force_coefficient_h)
        .def("calculate_force", &BasicNodePatch::calculate_force);

    class_<BasicTriangle>("BasicTriangle", init<const Vector&, const Vector&, const Vector&>())
        .def("normal", &BasicTriangle::get_planar_normal, return_value_policy<copy_const_reference>())
        .def("area", &BasicTriangle::get_planar_area)
        .def("centroid", &BasicTriangle::get_planar_centroid, return_value_policy<copy_const_reference>());
        
    class_<Triangle>("Triangle", init<const Triangle&>())
        .def("centre", &Triangle::get_planar_centroid, return_value_policy<copy_const_reference>())
        .def("normal", &Triangle::get_planar_normal, return_value_policy<copy_const_reference>())
        .def("planar_area", &Triangle::get_planar_area)
        .def("v1", &Triangle::get_v1, return_value_policy<copy_const_reference>())
        .def("v2", &Triangle::get_v2, return_value_policy<copy_const_reference>())
        .def("v3", &Triangle::get_v3, return_value_policy<copy_const_reference>())
        .def("v1_idx", &Triangle::get_v1_idx)
        .def("v2_idx", &Triangle::get_v2_idx)
        .def("v3_idx", &Triangle::get_v3_idx)
        .def("calculate_force", &Triangle::calculate_force);

    class_<HybridBezierTriangle>("HybridBezierTriangle", init<const Triangle&>())
        .def("calculate_force", &HybridBezierTriangle::calculate_force);
    class_<PNG1_Triangle>("PNG1_Triangle", init<const Triangle&>())
        .def("calculate_force", &PNG1_Triangle::calculate_force);

    class_<Vertex>("Vertex", init<const Vertex&>())
        .def("vertex", &Vertex::get_vertex, return_value_policy<reference_existing_object>())
        .def("normal", &Vertex::get_normal, return_value_policy<reference_existing_object>());

    class_<Mesh>("_Mesh", init<const std::string&, const std::string&, bool>())
        .def(init<const std::string&>())
        .def(init<const Mesh&>())
        .def(init<>())
        .def("init", &Mesh::init)
        .def("get_triangle", &Mesh::get_triangle, return_value_policy<reference_existing_object>())
        .def("get_node_patch", &Mesh::get_node_patch, return_value_policy<reference_existing_object>())
        .def("get_charge", &Mesh::get_charge, return_value_policy<reference_existing_object>())
        .def("get_vertex", &Mesh::get_vertex, return_value_policy<reference_existing_object>())
        .def("add_vertex", &Mesh::add_vertex)
        .def("define_triangle", &Mesh::define_triangle)
        .def_readonly("num_triangles", &Mesh::get_num_triangles)
        .def_readonly("num_node_patches", &Mesh::get_num_node_patches)
        .def_readonly("num_vertices", &Mesh::get_num_node_patches) // number of vertices sane as num node patches
        .def_readonly("num_charges", &Mesh::get_num_charges) 
        .def("get_centre", &Mesh::get_centre, return_value_policy<copy_const_reference>())
        .def("set_centre", &Mesh::set_centre)
        .def("init_energy_precalcs", &Mesh::init_energy_precalcs)
        .def("get_radius", &Mesh::get_radius)
        .def("calculate_volume", &Mesh::calculate_volume, "Calculate the volume enclosed by the mesh.")
        .def("kinemage_fh_vals", &Mesh::kinemage_fh_vals)
        .def("calculate_forces", &Mesh::py_calculate_forces)
        .def("calculate_energy", &Mesh::py_calculate_energy)
        .def("kinemage_node_patches", &Mesh::kinemage_node_patches)
        .def("set_quad_points_per_triangle", &Mesh::set_quad_points_per_triangle)
        .def("set_qual_points_per_triangle", &Mesh::set_qual_points_per_triangle);
;

    class_< BasicOctree<Vector> >("Octree", init<unsigned int, const Vector&, double>())
        .def("insert", &BasicOctree<Vector>::insert)
        .def("nearest", &BasicOctree<Vector>::get_nearest, return_value_policy<copy_const_reference>())
        .def("finalize", &BasicOctree<Vector>::finalize)
        .def("spew", &BasicOctree<Vector>::spew)
        .def("rejig", &BasicOctree<Vector>::rejig)
        .def("check_for_duplicate", &BasicOctree<Vector>::check_for_duplicate)
        .def_readonly("num_octree_nodes", &BasicOctree<Vector>::num_octree_nodes)
        .def_readonly("num_explicit_calcs", &BasicOctree<Vector>::num_explicit_calcs)
        .def_readonly("max_octree_depth", &BasicOctree<Vector>::max_octree_depth)
        .def_readonly("num_fmm_evals", &BasicOctree<Vector>::num_fmm_evals)
        .def_readonly("num_interaction_list_moves", &BasicOctree<Vector>::num_interaction_list_moves);
        
    class_< BasicOctree<BasicNodePatch> >("NodePatchTree", init<unsigned int, const Vector&, double>())
        .def("insert", &BasicOctree<BasicNodePatch>::insert)
        .def("rejig", &BasicOctree<BasicNodePatch>::rejig)
        .def("nearest", &BasicOctree<BasicNodePatch>::get_nearest, return_value_policy<copy_const_reference>());

    class_<RunInfo>("RunInfo", init<>())
        .def("__str__", &RunInfo::str);
        
    class_<BEEP>("_BEEP", init<double, double, int, int, int, bool>())
        .def(init<const std::string&, bool, bool>())
        .def("load_library_mesh", &BEEP::py_load_library_mesh, return_value_policy<copy_non_const_reference>())
        .def("clear_mesh_instances", &BEEP::py_clear_mesh_instances)
        .def("insert_mesh_instance", &BEEP::py_insert_mesh_instance, return_value_policy<copy_non_const_reference>())
        .def("py_kinemage" , &BEEP::py_kinemage)
        .def("solve", &BEEP::solve)
        .def("benchmark", &BEEP::benchmark)
        .def("calculate_energies", &BEEP::calculate_energies)
        .def("calculate_forces", &BEEP::calculate_forces)
        .def("reset_library_fh_vals", &BEEP::reset_library_fh_vals)
        .def("write_fh", &BEEP::write_fh)
        .def_readonly("num_node_patches", &BEEP::get_total_patches)
        .def_readonly("bem_neighbourhood_size", &BEEP::get_bem_neighbourhood_size)
        .def("set_bem_neighbourhood_size", &BEEP::set_bem_neighbourhood_size)
        .def("get_patch", &BEEP::get_patch, return_value_policy<copy_const_reference>())
        .def_readonly("num_mesh_instances", &BEEP::get_num_mesh_instances)
        .def("get_mesh_instance", &BEEP::get_mesh_instance, return_value_policy<copy_non_const_reference>())
        .def("set_bounding_cube", &BEEP::set_bounding_cube)
        .def("get_potential_at_point", &BEEP::get_potential_at_point)
        .def("write_opendx_pts_xyz", &BEEP::write_opendx_xyz)
        .def("write_opendx", &BEEP::write_opendx)
        .def("write_matrix", &BEEP::write_matrix)
        .def("skip_precalcs", &BEEP::hack_skip_precalcs)
        .def("set_quad_points_per_triangle", &BEEP::set_quad_points_per_triangle)
        .def("set_qual_points_per_triangle", &BEEP::set_qual_points_per_triangle);
        
    class_<MeshInstance>("_MeshInstance", init<const MeshInstance&>())
        .def_readonly("num_node_patches", &MeshInstance::get_num_node_patches)
        .def_readonly("num_charges", &MeshInstance::get_num_charges)
        .def("get_node_patch", &MeshInstance::get_node_patch, return_value_policy<copy_non_const_reference>())
        .def("get_charge", &MeshInstance::get_charge, return_value_policy<copy_const_reference>())
        .def("get_mesh", &MeshInstance::get_ref_mesh, return_value_policy<copy_const_reference>())
        .def("set_dielectrics", &MeshInstance::set_dielectrics)
        .def("isSilent", &MeshInstance::isSilent)
        .def("set_silent", &MeshInstance::set_silent)
        .def("unset_silent", &MeshInstance::unset_silent)
        .def("get_local_centre", &MeshInstance::get_local_centre, return_value_policy<copy_const_reference>())
        .def("get_xyz_offset", &MeshInstance::get_xyz_offset, return_value_policy<copy_const_reference>())
        .def("get_rotation", &MeshInstance::get_rotation, return_value_policy<copy_const_reference>())
        .def("get_radius", &MeshInstance::get_radius)
        .def("get_ecm_charges", &MeshInstance::get_ecm_charges)
        .def("set_quad_points_per_triangle", &MeshInstance::set_quad_points_per_triangle)
        .def("set_qual_points_per_triangle", &MeshInstance::set_qual_points_per_triangle);
        
    class_<GridParms>("GridParms")
        .def(init<>())
        .def_readwrite("x_pts", &GridParms::x_pts)
        .def_readwrite("y_pts", &GridParms::y_pts)
        .def_readwrite("z_pts", &GridParms::z_pts)
        .def_readwrite("x_angstroms", &GridParms::x_angstroms)
        .def_readwrite("y_angstroms", &GridParms::y_angstroms)
        .def_readwrite("z_angstroms", &GridParms::z_angstroms)
        .def_readwrite("origin_x", &GridParms::origin_x)
        .def_readwrite("origin_y", &GridParms::origin_y)
        .def_readwrite("origin_z", &GridParms::origin_z)
        .def("deltax",&GridParms::deltax, "Get x increment")
        .def("deltay",&GridParms::deltay, "Get y increment")
        .def("deltaz",&GridParms::deltaz, "Get z increment")
        .def("num_grid_points",&GridParms::num_grid_points, "Return total number of grid points")
        .def("GridIdxToPos",&GridParms::GridIdxToPos, "Convert xyz grid indices to VectorPoint")
        .def("OpenDX_Suffix", &GridParms::OpenDX_Suffix, "Return the OpenDX suffix string")
        .def("OpenDX_Preamble", &GridParms::OpenDX_Preamble, "Return the OpenDX preamble");
}

