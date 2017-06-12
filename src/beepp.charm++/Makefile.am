include $(top_srcdir)/common.am
CHARMC = charmc 
CXX = $(CHARMC) -O2 -tracemode projections
#CXX = $(CHARMC)
AUTOMAKE_OPTIONS = no-dependencies
EXTRA_DIST = fh_values_nodegroup.ci gmres.ci iteration_handler.ci main.ci mesh_library.ci mesh_working_copy.ci opencl_nodegroup.ci rhs_handler.ci vanilla_fmm.ci bem_fmm.ci 
bin_PROGRAMS = beepp
nodist_beepp_SOURCES = main.decl.h \
		       main.def.h \
		       mesh_working_copy.def.h \
		       mesh_working_copy.decl.h \
		       mesh_library.def.h \
		       mesh_library.decl.h \
		       rhs_handler.def.h \
		       rhs_handler.decl.h \
		       gmres.def.h \
		       gmres.decl.h \
		       iteration_handler.def.h \
		       iteration_handler.decl.h \
		       fh_values_nodegroup.def.h \
		       fh_values_nodegroup.decl.h \
		       opencl_nodegroup.def.h \
		       opencl_nodegroup.decl.h \
			   vanilla_fmm.def.h \
			   vanilla_fmm.decl.h \
			   bem_fmm.def.h \
			   bem_fmm.decl.h

beepp_SOURCES = prerequisites.h \
				vanilla_fmm_prerequisites.h \
		        main.cpp \
                main.h \
                mesh_working_copy.cpp \
                mesh_working_copy.h \
                mesh_library.cpp \
                mesh_library.h \
                rhs_handler.cpp \
                rhs_handler.h \
                gmres.cpp \
                gmres.h \
                iteration_handler.cpp \
                iteration_handler.h \
                fh_values_nodegroup.cpp \
                fh_values_nodegroup.h \
                opencl_nodegroup.cpp \
                opencl_nodegroup.h \
				charm_node_patch.h \
				charm_node_patch.cpp \
				vanilla_fmm.h \
				vanilla_fmm.cpp \
				bem_fmm.h \
				bem_fmm.cpp
                
beepp_CPPFLAGS = -fopenmp
beepp_LDADD = ../bem/libbem_charmc.la ../external/tinyxml/libtinyxml.la ../fmm/libfmm.la ../opencl/libopencl_handler.la
beepp_LDFLAGS = -lgsl -lgslcblas -larchive -lboost_filesystem -lboost_thread -pg -lgomp -lboost_program_options

# Charm headers generated from .ci interface files
fh_values_nodegroup.def.h:
	${CHARMC} fh_values_nodegroup.ci
fh_values_nodegroup.decl.h:
	${CHARMC} fh_values_nodegroup.ci
fmm_worker.def.h:
	${CHARMC} fmm_worker.ci
fmm_worker.decl.h:
	${CHARMC} fmm_worker.ci
gmres.def.h:
	${CHARMC} gmres.ci
gmres.decl.h:
	${CHARMC} gmres.ci
iteration_handler.def.h:
	${CHARMC} iteration_handler.ci
iteration_handler.decl.h:
	${CHARMC} iteration_handler.ci
main.def.h:
	${CHARMC} main.ci
main.decl.h:
	${CHARMC} main.ci
mesh_library.def.h:
	${CHARMC} mesh_library.ci
mesh_library.decl.h:
	${CHARMC} mesh_library.ci
mesh_working_copy.def.h:
	${CHARMC} mesh_working_copy.ci
mesh_working_copy.decl.h:
	${CHARMC} mesh_working_copy.ci
opencl_nodegroup.def.h:
	${CHARMC} opencl_nodegroup.ci
opencl_nodegroup.decl.h:
	${CHARMC} opencl_nodegroup.ci
parallel_fmm_octree.def.h:
	${CHARMC} parallel_fmm_octree.ci
parallel_fmm_octree.decl.h:
	${CHARMC} parallel_fmm_octree.ci
rhs_handler.def.h:
	${CHARMC} rhs_handler.ci
rhs_handler.decl.h:
	${CHARMC} rhs_handler.ci
