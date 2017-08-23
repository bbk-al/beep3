/*      Author: david fallaize    Created on: 26 Jul 2010 */

/*! \file mesh_tarball.h
 * \brief This module declares the mesh tarball class providing utilities
 * for reading compressed tar files containing mesh files.
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
#include <boost/filesystem/fstream.hpp>  // Not needed by 1.62.0
#endif

namespace fs = boost::filesystem;	// Easier to swap to std::filesystem in '17

class MeshTarball {

public:

    MeshTarball() : mtz_filename("") {}
	MeshTarball(const std::string& filename)
    : mtz_filename{filename},
	  work_dir {fs::temp_directory_path().string()
				+ work_dir.preferred_separator
				+ fs::unique_path("beep-%%%%%%%%%%").string()}
	{
		init();
	}

    virtual ~MeshTarball();

    class MeshTarball_Exception : public std::exception
    {
    public:
        MeshTarball_Exception() : std::exception() {}

    };

    // these are mandatory elements
    fs::path get_xyzqr_filename() const {
		return prepend_work_dir(xyzqr_filename);
	}
    fs::path get_centre_filename() const {
		return prepend_work_dir(centre_filename);
	}
    fs::path get_energies_filename() const {
		return prepend_work_dir(energies_filename);
	}
    fs::path get_mesh_filename() const {
		return prepend_work_dir(mesh_filename);
	}
	fs::path get_fh_filename() const { return prepend_work_dir(fh_filename); }
	fs::path get_ecm_filename() const { return prepend_work_dir(ecm_filename); }
#ifndef PREVOLHE
	// Optional elements
    fs::path get_mesh2_filename() const {
		return get_opt_filename(mesh2_filename);
	}
    fs::path get_energies2_filename() const {
		return get_opt_filename(energies2_filename);
	}
#endif // PREVOLHE

private:
#ifndef PREVOLHE
    fs::path get_opt_filename(const fs::path& fn) const {
		if (fn.empty()) return fn;
		return prepend_work_dir(fn.string());
	}
#endif // PREVOLHE

    inline fs::path prepend_work_dir(const std::string& filename) const;

    void init();
    void extract(const std::string& filename);
    int  copy_data(struct archive *ar, struct archive *aw);
    void parse_definition_xml();

	// Attributes
    std::string mtz_filename;
	fs::path work_dir;

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
#ifndef PREVOLHE
    static const std::string MESH2_TAG;
    static const std::string ENERGIES2_TAG;
#endif // PREVOLHE

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
#ifndef PREVOLHE
    std::string mesh2_filename;
    std::string energies2_filename;
#endif // PREVOLHE
};

// inlined methods
inline fs::path
MeshTarball::prepend_work_dir(const std::string& filename) const
{
	if (filename.empty()) { throw MeshTarball_Exception(); }
	fs::path rv{work_dir};
	rv /= fs::path{filename};
	return rv;
}

#endif /* MESH_TARBALL_H_ */
