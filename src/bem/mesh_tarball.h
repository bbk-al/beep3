/*
* mesh_tarball.h
*
*  Created on: 26 Jul 2010
*      Author: david
*/

#ifndef MESH_TARBALL_H_
#define MESH_TARBALL_H_

#include <string>
#include <exception>

// boost filesystem for messing with directories
#include <boost/version.hpp>
#if BOOST_VERSION < 103400
// legion has practically the oldest version of boost in existence
// (ok not really, but it does need updating...)
#include "boost/filesystem/operations.hpp" 
#else
#define BOOST_FILESYSTEM_VERSION 3
#include <boost/filesystem.hpp>
#endif

class mesh_tarball {

public:


    mesh_tarball() : mtz_filename("") {}
    mesh_tarball(const std::string& filename);
    virtual ~mesh_tarball();

    void init(const std::string& filename);

    class Mesh_TarBall_Exception : public std::exception
    {
    public:
        Mesh_TarBall_Exception() : std::exception() {}

    };

    // these are mandatory elements
    inline std::string prepend_work_dir(const std::string& filename) const {
        if (filename.empty()) { throw Mesh_TarBall_Exception(); }
        return work_dir + "/" + filename;
    }
    inline std::string get_xyzqr_filename() const { return prepend_work_dir(xyzqr_filename); }
    inline std::string get_mesh_filename() const { return prepend_work_dir(mesh_filename); }
    inline std::string get_centre_filename() const { return prepend_work_dir(centre_filename); }
    inline std::string get_energies_filename() const { return prepend_work_dir(energies_filename); }
    inline std::string get_fh_filename() const 
    { 
        if (fh_filename.empty()) {
            return ""; 
        }
        return prepend_work_dir(fh_filename); 
    };
    
    inline std::string get_ecm_filename() const
    {
        if (ecm_filename.empty()) {
            return ""; 
        }
        return prepend_work_dir(ecm_filename); 

    };

private:

    std::string mtz_filename;
    std::string work_dir;
    void extract(const std::string& filename);
    int  copy_data(struct archive *ar, struct archive *aw);
    void parse_definition_xml();

    // compulsory elements
    static const std::string MESH_DEFINITION_XML;
    static const std::string MESH_DEFINITION_XML_ROOT_ELEMENT;
    static const std::string PQR_TAG;
    static const std::string MESH_TAG;
    static const std::string CENTRE_TAG;
    static const std::string ENERGIES_TAG;
    static const std::string XYZQR_TAG;

    // optional elements
    static const std::string PDB_TAG; // unused
    static const std::string VERTEX_NORMALS_TAG; // unused
    static const std::string FH_TAG;
    static const std::string TENSOR_TAG; // unused
    static const std::string ELLIPSOID_TAG; // unused
    static const std::string ECM_TAG;

    // internal storage of the filenames we find
    std::string pdb_filename;
    std::string pqr_filename;
    std::string xyzqr_filename;
    std::string vertex_normal_filename;
    std::string mesh_filename;
    std::string centre_filename;
    std::string energies_filename;
    std::string fh_filename;
    std::string diff_tensor_filename;
    std::string ellipsoid_filename;
    std::string ecm_filename;

};

#endif /* MESH_TARBALL_H_ */
