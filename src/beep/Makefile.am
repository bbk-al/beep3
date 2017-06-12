include $(top_srcdir)/common.am
bin_PROGRAMS = beep gb
beep_SOURCES = beep.cpp beep.h beep_main.cpp
beep_LDADD = ../bem/libbem.la 
beep_LDFLAGS = -pg -lgomp -lboost_program_options
beep_CPPFLAGS = -I/usr/local/cuda/include -g -fopenmp

gb_SOURCES = generalized_born.cpp generalized_born.h ../bem/triangle.cpp ../bem/vertex.cpp 
gb_LDFLAGS = -lgsl -lgslcblas -lboost_filesystem -lboost_system 
gb_CPPFLAGS = -g
