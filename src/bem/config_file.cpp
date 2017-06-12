/*
* config_file.cpp
*
*  Created on: 23 Jul 2010
*      Author: david
*/

#include <cassert>
#include <iostream>
#include "config_file.h"
#include <sstream>
#include <boost/lexical_cast.hpp>

// boost filesystem for messing with directories
#include <boost/version.hpp>
#if BOOST_VERSION < 103400
// legion has practically the oldest version of boost in existence
// (ok not really, but it does need updating...)
#include "boost/filesystem/path.hpp" 
#else
#define BOOST_FILESYSTEM_VERSION 3
#include <boost/filesystem.hpp>
#endif

// Use the TinyXML library for XML parsing -- see external folder for details
#include "../external/tinyxml/tinyxml.h"

const std::string ConfigFile::XML_CONFIG_ROOT_ELEMENT = "BEM_Config";

void ConfigFile::init(const std::string& filename)
{
    
#if BOOST_VERSION < 103400
    std::string root_path = boost::filesystem::path(filename).branch_path().string() + "/";
#else
    std::string root_path = boost::filesystem::path(filename).parent_path().string() + "/";
#endif
    
    
    // Open the xml config file and parse it
    // Populate the ConfigFile members from the xml

    TiXmlDocument doc( filename );
    bool loadOkay = doc.LoadFile();

    if (!loadOkay) {
        std::cout << "Failed to load file: " << filename << std::endl;
        throw ParseException();
    }

    // Loaded config file, now lets parse the elements
    TiXmlHandle hDoc(&doc);
    TiXmlHandle hRoot(0);

    // block: name
    {
        TiXmlElement* pElem = hDoc.FirstChildElement().Element();
        // should always have a valid root but handle gracefully if it does
        if (!pElem) return;
        const std::string bem_config_name = pElem->ValueStr();
        if (bem_config_name != XML_CONFIG_ROOT_ELEMENT)
        {
            std::cerr << "Root element should be: " << XML_CONFIG_ROOT_ELEMENT << " (got: " << bem_config_name << ")" << std::endl;
            throw ParseException();
        }

        // save this for later
        hRoot=TiXmlHandle(pElem);

        // Global Quad/Qual (can be overridden within each instance elements)
        const char *pGlobalQuad = pElem->Attribute("quad");
        const char *pGlobalQual = pElem->Attribute("qual");
        global_default_quads = (pGlobalQuad) ? boost::lexical_cast<int>(std::string(pGlobalQuad)) : -1;
        global_default_quals = (pGlobalQual) ? boost::lexical_cast<int>(std::string(pGlobalQual)) : -1;

    }

    {
        TiXmlElement* pElem = hRoot.FirstChild( "solvent" ).Element();
        if (pElem == NULL)
        {
            std::cerr << "ERROR: Failed to find <solvent> tag (in which solvent dielectric & kappa are required)." << std::endl;
            throw ParseException();
        }

        const char *pDielectric = pElem->Attribute("dielectric");
        if (pDielectric == NULL)
        {
            std::cerr << "ERROR: Failed to find a value for solvent_dielectric in solvent definition." << std::endl;
            throw ParseException();
        }
        solvent_dielectric = boost::lexical_cast<double>(std::string(pDielectric));
        
        
        const char *pKappa = pElem->Attribute("kappa");
        solvent_kappa = boost::lexical_cast<double>(std::string(pKappa));
        if (pKappa == NULL)
        {
            std::cerr << "ERROR: Failed to find a value for kappa in solvent definition." << std::endl;
            throw ParseException();
        }


    }

    // output filename
    {
        TiXmlElement* pElem = hRoot.FirstChild( "output" ).Element();
        if (pElem == NULL) {
            std::cerr << "Failed to find an output tag in " << filename << std::endl;
            throw ParseException();
        }
        const char *pText=pElem->GetText();
        output_file = pText;
    }

    // block: Mesh library
    mesh_library.clear();
    {
        TiXmlElement* pElem = hRoot.FirstChild( "library" ).FirstChild().Element();
        for( ; pElem; pElem=pElem->NextSiblingElement())
        {
            const char *pKey=pElem->Value();
            const char *pText=pElem->GetText();
            const char *pId = pElem->Attribute("id");
            assert(pKey && pText && pId);
            std::string mesh_str(pKey);
            if (mesh_str != "mesh")
            {
                std::cerr << "Expected <mesh> tags within <library> block." << std::endl;
                throw ParseException();
            }

            unsigned int id = boost::lexical_cast<unsigned int>(std::string(pId));

            // check that mesh id is consistent with ordering within the xml
            if (id != mesh_library.size())
            {
                std::cerr << "Mesh id's must correspond to ordering within <library> block," << std::endl;
                throw ParseException();
            }
            mesh_library.push_back(CFG_MeshDescription(root_path + std::string(pText), id));

        }
    }

    layout.clear();
    {
        TiXmlElement* pElem = hRoot.FirstChild( "layout" ).FirstChild().Element();
        assert(pElem != NULL);
        for( ; pElem; pElem=pElem->NextSiblingElement())
        {
            const char *pKey=pElem->Value();
            const char *pInstanceId = pElem->Attribute("instance_id");
            const char *pMeshId = pElem->Attribute("mesh_id");
            const char *pProteinDielectric = pElem->Attribute("dielectric");
            const char *pX = pElem->Attribute("x");
            const char *pY = pElem->Attribute("y");
            const char *pZ = pElem->Attribute("z");
            const char *pA = pElem->Attribute("a");
            const char *pB = pElem->Attribute("b");
            const char *pC = pElem->Attribute("c");
            const char *pD = pElem->Attribute("d");

            std::string mesh_str(pKey);
            if (mesh_str != "instance")
            {
                std::cerr << "Expected <instance> tags within <layout> block." << std::endl;
                throw ParseException();
            }

            if (!(pInstanceId && pMeshId))
            {
                std::cerr << "Missing id attributes in <instance> tag." << std::endl;
                throw ParseException();
            }
            unsigned int instance_id = boost::lexical_cast<unsigned int>(std::string(pInstanceId));
            unsigned int mesh_id = boost::lexical_cast<unsigned int>(std::string(pMeshId));
            double protein_dielectric;

            if (pProteinDielectric)
            {
                protein_dielectric = boost::lexical_cast<double>(std::string(pProteinDielectric));
            }
            else
            {
                std::cerr << "Missing dielectric attribute in <instance> tag." << std::endl;
                throw ParseException();
            }

            if (!(pX && pY && pZ))
            {
                std::cerr << "Missing x/y/z attributes in <instance> tag." << std::endl;
                throw ParseException();
            }

            double x = boost::lexical_cast<double>(std::string(pX));
            double y = boost::lexical_cast<double>(std::string(pY));
            double z = boost::lexical_cast<double>(std::string(pZ));

            double a,b,c,d;
            if (!(pA && pB && pC && pD))
            {
                // default zero rotation quaternion
                a = 1.0;
                b = 0.0;
                c = 0.0;
                d = 0.0; 
            }
            else
            {
                a = boost::lexical_cast<double>(std::string(pA));
                b = boost::lexical_cast<double>(std::string(pB));
                c = boost::lexical_cast<double>(std::string(pC));
                d = boost::lexical_cast<double>(std::string(pD));
                
                // hacky fix for zero quaternions- make the real part equal 1
                if (a==0 && b==0 && c==0 && d==0) { a = 1.0; }
            }

            Vector xyz_offset(x,y,z);
            Quaternion rotation(a,b,c,d);
            if (fabs(rotation.norm() - 1.0) > 1e-3) {
                std::cerr << "Quaternion is not close to being a unit quaternion --> not a valid rotation?" << std::endl;
                throw ParseException();
            }
            rotation.normalise(); // make unit quaternion as close as possible with double precision fp...

            if (instance_id != layout.size())
            {
                std::cerr << "Instance id's must correspond to ordering within <layout> block." << std::endl;
                throw ParseException();
            }

            if (mesh_id >= mesh_library.size())
            {
                std::cerr << "Invalid mesh_id in <instance>." << std::endl;
                throw ParseException();
            }
            
            // Individual quad/qual specification (overrides global settings)
            const char *pQuad = pElem->Attribute("quad");
            const char *pQual = pElem->Attribute("qual");
            int quad_pts = (pQuad) ? boost::lexical_cast<int>(std::string(pQuad)) : global_default_quads;
            int qual_pts = (pQual) ? boost::lexical_cast<int>(std::string(pQual)) : global_default_quals;
            
            // create mesh_layout within config
/*            std::cout << "Read instance in layout: " << instance_id << " "
                    << mesh_id << " " << xyz_offset << " " << rotation << std::endl;*/
            layout.push_back( mesh_layout(instance_id, mesh_id, protein_dielectric, xyz_offset, rotation, quad_pts, qual_pts) );

        }
    }

    return;

}

