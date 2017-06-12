include $(top_srcdir)/common.am
lib_LTLIBRARIES = libBEEP.la
libBEEP_la_SOURCES = pybeep.cpp ../beep/beep.cpp 
libBEEP_la_LIBADD = ../bem/libbem.la 
libBEEP_la_LDFLAGS = $(LDPYBEEP)
libBEEP_la_CPPFLAGS = -fopenmp
EXTRA_DIST = add_results_to_mtz.py \
beep-ecm.py \
beep-nocrowding.py \
beep.py \
pybeep.py \
born_ion.py \
calc.py \
calc_energy.py \
calc_forces.py \
constants.py \
geometry.py \
gts_utils.py \
kintools.py \
mesh_utils.py \
opendx_utils.py \
pdb.py \
pqr2xyzqr.py \
pqrtools.py \
prepare.py \
print_mesh_stats.py 

