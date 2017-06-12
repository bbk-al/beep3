/*
* config_file.h
*
*  Created on: 23 Jul 2010
*      Author: david
*/

#ifndef CONFIG_FILE_H_
#define CONFIG_FILE_H_

#include <string>
#include <vector>
#include "../common/math_vector.h"

class ConfigFile {

public:

    struct CFG_MeshDescription
    {
        CFG_MeshDescription(const std::string& mesh, unsigned int id) : mesh_filename(mesh), mesh_id(id) {}
        CFG_MeshDescription(const char* mesh, unsigned int id) : mesh_filename(mesh), mesh_id(id) {}
        std::string mesh_filename;
        unsigned int mesh_id;
    };

    class mesh_layout
    {
    public:
        mesh_layout(unsigned int _instance_id,
                    unsigned int _mesh_id,
                    double _dielectric,
                    const Vector& _offset,
                    const Quaternion& _rotation,
                    int _quad_points,
                    int _qual_points) :
                     instance_id(_instance_id),
                     mesh_id(_mesh_id),
                     dielectric(_dielectric),
                     offset(_offset),
                     rotation(_rotation),
                     quad_points(_quad_points),
                     qual_points(_qual_points) {}

        unsigned int instance_id;
        unsigned int mesh_id;
        double dielectric;
        Vector offset;
        Quaternion rotation;
        int quad_points;
        int qual_points;
    };

    typedef std::vector<CFG_MeshDescription> MeshLibrary;
    typedef std::vector< mesh_layout > Layout;

    ConfigFile() : global_default_quads(-1), global_default_quals(-1) {}
    ConfigFile(const std::string& filename)  : global_default_quads(-1), global_default_quals(-1) { init(filename); }
    void init(const std::string& filename);

    MeshLibrary mesh_library;
    Layout layout;

    double solvent_kappa;
    double solvent_dielectric;
    unsigned int target_bem_explicits;
    std::string output_file;
    int global_default_quads;
    int global_default_quals;


    class ParseException : public std::exception
    {
    public:
        ParseException() : std::exception() {}

    };

    static const std::string XML_CONFIG_ROOT_ELEMENT;

};



#endif /* CONFIG_FILE_H_ */
