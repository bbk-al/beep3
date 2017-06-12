include $(top_srcdir)/common.am
noinst_LTLIBRARIES = libopencl_handler.la
libopencl_handler_la_SOURCES = opencl_handler.cpp opencl_handler.h opencl_workblob.h
libopencl_handler_la_CPPFLAGS = @CL_CFLAGS@
libopencl_handler_la_LDFLAGS = -fPIC @CL_LIBS@ -lboost_thread
