/*      Author: david fallaize   Created on: 26 Jul 2010 */

/*! \file mesh_tarball.cpp
 * \brief This module implements the mesh tarball class.
 */

#include "mesh_tarball.h"

#include <sstream>
#include <iostream>
#include <exception>
#include <vector>

// libarchive extraction includes
#include <sys/types.h>
#include <sys/stat.h>
#include <archive.h>
#include <archive_entry.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>


// tinyxml for parsing the xml description file
#include "../external/tinyxml/tinyxml.h"
#include <boost/tuple/tuple.hpp> // handy!

namespace fs = boost::filesystem;	// Easier to swap to std::filesystem in '17

// static const's for xml tags
const std::string MeshTarball::MESH_DEFINITION_XML = "definition.xml";
const std::string MeshTarball::MESH_DEFINITION_XML_ROOT_ELEMENT = "contents";
const std::string MeshTarball::PDB_TAG = "pdb";
const std::string MeshTarball::PQR_TAG = "pqr";
const std::string MeshTarball::XYZQR_TAG = "xyzqr";
const std::string MeshTarball::VERTEX_NORMALS_TAG = "vertex_normals";
const std::string MeshTarball::MESH_TAG = "mesh";
const std::string MeshTarball::ENERGIES_TAG = "energies";
const std::string MeshTarball::FH_TAG = "fh";
const std::string MeshTarball::TENSOR_TAG = "diffusion_tensor";
const std::string MeshTarball::CENTRE_TAG = "centre";
const std::string MeshTarball::ELLIPSOID_TAG = "ellipsoid";
const std::string MeshTarball::ECM_TAG = "ecm";


void MeshTarball::init() {
	// first change directory to a safe working location
#ifdef __DELETED__
    char work_dir_template[] = "tmpdir.XXXXXX";
    char* d = (mkdtemp(work_dir_template));
    if (d != nullptr) {
        work_dir = std::string(d);
    }
    else {
        std::cerr << "ERROR: Failed to create temporary working folder."
		          << std::endl;
        throw MeshTarball_Exception();
    }
#else // __DELETED__
	fs::create_directories(work_dir);
#endif // __DELETED__

//     try 
//     {
//         work_dir = "tmpdir." + filename;
//         boost::filesystem::create_directory(work_dir);
//     }
//     catch (...) {
//         std::cerr << "ERROR: Failed to create temporary working folder." << std::endl;
//         throw MeshTarball_Exception();
//     }

    // extract the contents of the mesh tarball
#ifdef __DELETED__
    chdir(work_dir.c_str());
    extract("../" + mtz_filename);
    chdir("..");
#else // __DELETED__
	fs::path cwd{fs::current_path()};
	fs::path mtz{mtz_filename};
	if (mtz.is_relative()) mtz = cwd / mtz;
	fs::current_path(work_dir);
	extract(mtz.string());
	fs::current_path(cwd);
#endif // __DELETED__

    // try to parse "definition.xml"
    parse_definition_xml();
}

MeshTarball::~MeshTarball() {
	
#ifndef __DELETED__
	// Directories can be left behind because of temporary nfs files
	// which disappear after a short delay, allowing the top level
	// directory to be deleted only after this process exits.
	do {  // Do-once loop to allow break out if something works
		try {
			// delete the temporary folder
			fs::remove_all(work_dir);
			break;
		}
		catch (std::exception& e) {
			std::cerr << "MeshTarball: unable to remove temporary directory "
			          << work_dir << ": " << e.what() << std::endl;
			// But note this error is fixed by using temp_directory_path above:
			std::cerr << "This may be due to a .nfs file which will "
			             "disappear when this process exits." << std::endl;
			std::cerr << "List of entries:" << std::endl;
			for (fs::directory_iterator itr(work_dir), end; itr != end; ++itr) {
				std::cerr << itr->path() << std::endl;
			}
			std::cerr << "End of entries" << std::endl;
		} // leaving the directory still in place
	} while (false);
#else // __DELETED__
    try {
        // delete the temporary folder
        for (fs::directory_iterator itr(work_dir), end; itr != end; ++itr) {
            fs::remove(*itr);
        }
        fs::remove(work_dir);
        
        
    }
    catch (...) {}
#endif // __DELETED__
}

int MeshTarball::copy_data(struct archive *ar, struct archive *aw)
{
    int r;
    const void *buff;
    size_t size;
    off_t offset;

    for (;;) {
        r = archive_read_data_block(ar, &buff, &size, &offset);
        if (r == ARCHIVE_EOF) return (ARCHIVE_OK);
        if (r != ARCHIVE_OK) return (r);

        r = archive_write_data_block(aw, buff, size, offset);
        if (r != ARCHIVE_OK) return (r);
    }
    return r;
}

void MeshTarball::extract(const std::string& filename)
{
    struct archive *a;
    struct archive *ext;
    struct archive_entry *entry;
    int r;

    a = archive_read_new();
    ext = archive_write_disk_new();
    archive_write_disk_set_options(ext, ARCHIVE_EXTRACT_SECURE_NODOTDOT);
    archive_read_support_filter_bzip2(a);
    archive_read_support_filter_gzip(a);
    archive_read_support_filter_compress(a);
    archive_read_support_format_tar(a);
    archive_read_support_format_cpio(a);
    //archive_write_disk_set_standard_lookup(ext);

    if ((r = archive_read_open_filename(a, filename.c_str(), 10240))) {
        std::cerr << archive_error_string(a) << std::endl;
        throw MeshTarball_Exception();
    }
    for (;;) {
        r = archive_read_next_header(a, &entry);
        if (r == ARCHIVE_EOF) break;

        if (r != ARCHIVE_OK) {
            std::cerr << archive_error_string(a) << std::endl;
            throw MeshTarball_Exception();
        }
        r = archive_write_header(ext, entry);
        if (r != ARCHIVE_OK) {
            std::cerr << archive_error_string(a) << std::endl;
            throw MeshTarball_Exception();
        }
        else {
            copy_data(a, ext);
        }
    }
    archive_read_close(a);
    archive_read_free(a);
}

void MeshTarball::parse_definition_xml()
{

    // Open the xml config file and parse it
#ifdef __DELETED__
    std::string xml = work_dir + "/" + MESH_DEFINITION_XML;
#else
    std::string xml = work_dir.string() + "/" + MESH_DEFINITION_XML;
#endif // __DELETED__
    TiXmlDocument doc(xml.c_str());
    bool loadOkay = doc.LoadFile();
    bool failed = false;

    if (!loadOkay) {
        std::cerr << "ERROR: Failed to load " << MESH_DEFINITION_XML
		          << " within mesh " << mtz_filename << std::endl;
        throw MeshTarball_Exception();
    }

    std::string err_fname_path = MESH_DEFINITION_XML + " ("
								+ mtz_filename + ")";

    // Loaded config file, now lets parse the elements
    TiXmlHandle hDoc(&doc);
    TiXmlElement* pElem;
    TiXmlHandle hRoot(0);

    // block: name
    {
        pElem=hDoc.FirstChildElement().Element();
        // should always have a valid root but handle gracefully if it does
        if (!pElem) {
            std::cerr << "ERROR: Unable to find root element in "
			          << err_fname_path << std::endl;
            throw MeshTarball_Exception();
        }
        std::string bem_config_name = pElem->Value();
        if (bem_config_name != MESH_DEFINITION_XML_ROOT_ELEMENT) {
            std::cerr << "ERROR: Root element should be: "
			          << MESH_DEFINITION_XML_ROOT_ELEMENT << " in "
			          << err_fname_path <<  std::endl;
            throw MeshTarball_Exception();
        }

        // save this for later
        hRoot=TiXmlHandle(pElem);
    }

    // a triplet is a tuple of 3 values - the third value (a bool) is whether
    // or not the tag is mandatory within the xml definition
    typedef boost::tuple<const std::string*, std::string*, bool> triplet;
    std::vector<triplet> key_val_defs;
    key_val_defs.push_back( triplet(&PDB_TAG, &pdb_filename, false) );
    key_val_defs.push_back( triplet(&PQR_TAG, &pqr_filename, false) );
    key_val_defs.push_back( triplet(&XYZQR_TAG, &xyzqr_filename, true) );
    key_val_defs.push_back( triplet(&VERTEX_NORMALS_TAG,
									&vertex_normal_filename, false) );
    key_val_defs.push_back( triplet(&MESH_TAG, &mesh_filename, true) );
    key_val_defs.push_back( triplet(&ENERGIES_TAG, &energies_filename, true) );
    key_val_defs.push_back( triplet(&MESH_TAG, &mesh_filename, true) );
    key_val_defs.push_back( triplet(&FH_TAG, &fh_filename, false) );
    key_val_defs.push_back( triplet(&TENSOR_TAG, &diff_tensor_filename, false)
						  );
    key_val_defs.push_back( triplet(&CENTRE_TAG, &centre_filename, true) );
    key_val_defs.push_back( triplet(&ELLIPSOID_TAG, &ellipsoid_filename, false)
						  );
    key_val_defs.push_back( triplet(&ECM_TAG, &ecm_filename, false) );

    // iterate over this list of allowable key/vals and fill in the struct
    for (std::vector<triplet>::iterator
			it=key_val_defs.begin(), end=key_val_defs.end();
		 it != end; ++it)
    {
        triplet& def = *it;
        const std::string& tag = *(def.get<0>());
        std::string& text = *(def.get<1>());
        bool mandatory = def.get<2>();

        pElem=hRoot.FirstChild( tag.c_str() ).Element();
        if (pElem) {
            text = std::string( pElem->GetText() );

            //if (pElem->FirstChild()) {
            //	std::cerr << "ERROR: Unexpected children for \"" << tag << "\" tag in " << err_fname_path << std::endl;
            //	failed = true;
            //}

            if(pElem->NextSiblingElement( tag.c_str() )) {
                std::cerr << "ERROR: Multiple \"" << tag << "\" tags in "
				          << err_fname_path << std::endl;
                failed = true;
            }

            // TODO: stat the file to check it really does exist in the meshtar
        }
        else {
            text = "";

            // If it is mandatory then failing to find the tag is a fatal
            // error -- cannot make a mesh without this information.
            // (e.g. gts/off file, charge locations)
            if (mandatory) {
                std::cerr << "ERROR: Failed to find mandatory tag \"" << tag
				          << "\" in " << err_fname_path << std::endl;
                failed = true;
            }
        }
    }

    if (failed) {
        throw MeshTarball_Exception();
    }
}

