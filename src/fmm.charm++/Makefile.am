include $(top_srcdir)/common.am
CHARMC = charmc 
CXX = $(CHARMC) -O2 -tracemode projections 
AUTOMAKE_OPTIONS = no-dependencies
EXTRA_DIST = main.ci opencl_nodegroup.ci vanilla_fmm.ci 
bin_PROGRAMS = pfmm
nodist_pfmm_SOURCES = main.decl.h \
		       main.def.h \
         	   vanilla_fmm.def.h \
		       vanilla_fmm.decl.h \
		       opencl_nodegroup.def.h \
		       opencl_nodegroup.decl.h

pfmm_SOURCES = vanilla_fmm_prerequisites.h \
  		       main.cpp \
               main.h \
               vanilla_fmm.cpp \
               vanilla_fmm.h \
               opencl_nodegroup.cpp \
               opencl_nodegroup.h
                
pfmm_LDADD = ../fmm/libfmm.la ../opencl/libopencl_handler.la
pfmm_LDFLAGS = -lgomp
pfmm_CPPFLAGS = -fopenmp 

