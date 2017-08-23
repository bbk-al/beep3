include $(top_srcdir)/common.am

noinst_LTLIBRARIES = libbem.la 
if WITH_CHARM_PLUS_PLUS 
noinst_LTLIBRARIES += libbem_charmc.la
endif

libbem_charmc_la_SOURCES = bem_kernels.cpp \
		config_file.cpp \
		mesh.cpp \
		meshing.cpp \
		mesh_instance.cpp \
		mesh_tarball.cpp \
		node_patch.cpp \
		spharm.cpp \
		triangle.cpp \
		vertex.cpp \
		gts_utils.h \
		off_utils.h \
		bem_kernels.h \
		bezier.h \
		png1.h \
		config_file.h \
		constants.h \
		edge.h \
		local_integrations.h \
		mesh.h \
		meshing.h \
		mesh_instance.h \
		mesh_tarball.h \
		node_patch.h \
		quad_point.h \
		spharm.h \
		triangle.h \
		vertex.h \
		opencl_bem.h \
		opencl_bem.cpp \
		opencl_bem_structs.h \
		opencl_bem_kernels.cl.h \
		opendx_utils.h

libbem_charmc_la_CPPFLAGS = -D__CHARMC__  
libbem_charmc_la_LIBADD = ../fmm/libfmm.la ../external/tinyxml/libtinyxml.la -lgsl -lgslcblas -lboost_filesystem -larchive -lboost_thread -lboost_system

libbem_la_SOURCES = bem_kernels.cpp \
		config_file.cpp \
		mesh.cpp \
		meshing.cpp \
		mesh_instance.cpp \
		mesh_tarball.cpp \
		node_patch.cpp \
		spharm.cpp \
		triangle.cpp \
		vertex.cpp \
		gts_utils.h \
		off_utils.h \
		bem_kernels.h \
		bezier.h \
		png1.h \
		config_file.h \
		constants.h \
		edge.h \
		local_integrations.h \
		mesh.h \
		meshing.h \
		mesh_instance.h \
		mesh_tarball.h \
		node_patch.h \
		quad_point.h \
		spharm.h \
		triangle.h \
		vertex.h \
		opencl_bem.h \
		opencl_bem.cpp \
		opencl_bem_structs.h \
		opencl_bem_kernels.cl.h \
		opendx_utils.h

libbem_la_LIBADD = ../fmm/libfmm.la ../external/tinyxml/libtinyxml.la ../opencl/libopencl_handler.la 
libbem_la_LDFLAGS = -lgsl -lgslcblas -lboost_filesystem -larchive -lboost_thread -lboost_system 
EXTRA_DIST = opencl_bem_kernels.cl

