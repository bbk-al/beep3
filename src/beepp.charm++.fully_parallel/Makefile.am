include $(top_srcdir)/common.am
CHARMC = charmc 
CXX = $(CHARMC) -O2 -tracemode projections -modules comlib
#CXX = $(CHARMC)
AUTOMAKE_OPTIONS = no-dependencies
EXTRA_DIST = explicit_worker.ci fh_values_nodegroup.ci fmm_globals_nodegroup.ci fmm_worker.ci gmres.ci iteration_handler.ci main.ci mesh_library.ci mesh_working_copy.ci opencl_nodegroup.ci rhs_handler.ci parallel_fmm_octree.ci vanilla_fmm_worker.ci vanilla_fmm_evals.ci
bin_PROGRAMS = beepp
nodist_beepp_SOURCES = main.decl.h \
		       main.def.h \
		       mesh_working_copy.def.h \
		       mesh_working_copy.decl.h \
		       fmm_globals_nodegroup.def.h \
		       fmm_globals_nodegroup.decl.h \
		       mesh_library.def.h \
		       mesh_library.decl.h \
		       parallel_fmm_octree.def.h \
		       parallel_fmm_octree.decl.h \
		       fmm_worker.def.h \
		       fmm_worker.decl.h \
		       rhs_handler.def.h \
		       rhs_handler.decl.h \
		       gmres.def.h \
		       gmres.decl.h \
		       iteration_handler.def.h \
		       iteration_handler.decl.h \
		       explicit_worker.def.h \
		       explicit_worker.decl.h \
		       fh_values_nodegroup.def.h \
		       fh_values_nodegroup.decl.h \
		       opencl_nodegroup.def.h \
		       opencl_nodegroup.decl.h \
			   vanilla_fmm_evals.def.h \
			   vanilla_fmm_evals.decl.h \
			   vanilla_fmm_worker.def.h \
			   vanilla_fmm_worker.decl.h

beepp_SOURCES = prerequisites.h \
				vanilla_fmm_prerequisites.h \
		        main.cpp \
                main.h \
                mesh_working_copy.cpp \
                mesh_working_copy.h \
                fmm_globals_nodegroup.cpp \
                fmm_globals_nodegroup.h \
                mesh_library.cpp \
                mesh_library.h \
                parallel_fmm_octree.cpp \
                parallel_fmm_octree.h \
                fmm_worker.cpp \
                fmm_worker.h \
                rhs_handler.cpp \
                rhs_handler.h \
                gmres.cpp \
                gmres.h \
                iteration_handler.cpp \
                iteration_handler.h \
                explicit_worker.cpp \
                explicit_worker.h \
                fh_values_nodegroup.cpp \
                fh_values_nodegroup.h \
                opencl_nodegroup.cpp \
                opencl_nodegroup.h \
				vanilla_fmm_evals.cpp \
				vanilla_fmm_evals.h \
				vanilla_fmm_worker.cpp \
				vanilla_fmm_worker.h \
				charm_node_patch.h \
				charm_node_patch.cpp
                
beepp_CPPFLAGS =
beepp_LDADD = ../bem/libbem_charmc.la ../external/tinyxml/libtinyxml.la ../fmm/libfmm.la ../opencl/libopencl_handler.la
beepp_LDFLAGS = -lgsl -lgslcblas -larchive -lboost_filesystem -lboost_thread -pg  

# Charm headers generated from .ci interface files
explicit_worker.def.h:
	${CHARMC} explicit_worker.ci
explicit_worker.decl.h:
	${CHARMC} explicit_worker.ci
fh_values_nodegroup.def.h:
	${CHARMC} fh_values_nodegroup.ci
fh_values_nodegroup.decl.h:
	${CHARMC} fh_values_nodegroup.ci
fmm_globals_nodegroup.def.h:
	${CHARMC} fmm_globals_nodegroup.ci
fmm_globals_nodegroup.decl.h:
	${CHARMC} fmm_globals_nodegroup.ci
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
