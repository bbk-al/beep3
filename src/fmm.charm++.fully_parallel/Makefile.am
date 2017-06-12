include $(top_srcdir)/common.am
CHARMC = charmc 
CXX = $(CHARMC) -O2 -tracemode projections -balancer RefineLB -modules commlib
AUTOMAKE_OPTIONS = no-dependencies
EXTRA_DIST = fmm_globals_nodegroup.ci vanilla_fmm_worker.ci main.ci opencl_nodegroup.ci parallel_fmm_octree.ci vanilla_fmm_evals.ci
bin_PROGRAMS = pfmm
nodist_pfmm_SOURCES = main.decl.h \
		       main.def.h \
		       fmm_globals_nodegroup.def.h \
		       fmm_globals_nodegroup.decl.h \
  	       	   parallel_fmm_octree.def.h \
		       parallel_fmm_octree.decl.h \
		       vanilla_fmm_worker.def.h \
		       vanilla_fmm_worker.decl.h \
			   vanilla_fmm_evals.def.h \
			   vanilla_fmm_evals.decl.h \
		       opencl_nodegroup.def.h \
		       opencl_nodegroup.decl.h

pfmm_SOURCES = vanilla_fmm_prerequisites.h \
  		       main.cpp \
               main.h \
               fmm_globals_nodegroup.cpp \
               fmm_globals_nodegroup.h \
               parallel_fmm_octree.cpp \
               parallel_fmm_octree.h \
               vanilla_fmm_worker.cpp \
               vanilla_fmm_worker.h \
			   vanilla_fmm_evals.cpp \
			   vanilla_fmm_evals.h \
               opencl_nodegroup.cpp \
               opencl_nodegroup.h
                
pfmm_LDADD = ../fmm/libfmm.la ../opencl/libopencl_handler.la
pfmm_LDFLAGS = -lgsl -lgslcblas -larchive -lboost_filesystem -lboost_thread  

