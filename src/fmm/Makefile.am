include $(top_srcdir)/common.am
noinst_LTLIBRARIES = libfmm.la
libfmm_la_SOURCES = fmm_bem_hybrid.h \
				    fmm.h \
					fmm_globals.h \
					fmm_octree.h \
					octree.h \
					fmm_math_funcs.h \
					bessel_gamma.h \
					eval_pt.h \
					opencl_fmm.h \
					opencl_fmm.cpp \
					opencl_fmm_kernels.cl.h \
					interaction_list.h \
					eval_pt.h \
					fmm_time_info.h
EXTRA_DIST = opencl_fmm_kernels.cl fmm.test gentest.py
libfmm_la_CPPFLAGS = -g -O2 -fopenmp
libfmm_la_LDFLAGS = -lgomp

bin_PROGRAMS = test_fmm
test_fmm_SOURCES = test_fmm.cpp eval_pt.h
test_fmm_LDADD = libfmm.la ../opencl/libopencl_handler.la
test_fmm_LDFLAGS = -pg -lgomp
test_fmm_CPPFLAGS = -g -O2 -fopenmp


