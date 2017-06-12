AC_INIT([beep], [0.2], [d.fallaize@ucl.ac.uk])
LT_INIT
AM_INIT_AUTOMAKE([-Wall -Werror])
AC_PROG_CXX
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile \
		 src/common/Makefile \
		 src/opencl/Makefile \
		 src/bem/Makefile \
		 src/fmm/Makefile \
		 src/fmm.charm++/Makefile \
		 src/beep/Makefile \
         src/beepp.charm++/Makefile \
		 src/pybeep/Makefile \
		 src/external/Makefile \
		 src/external/tinyxml/Makefile])
#AC_CONFIG_SUBDIRS([src/external/ltl-1.9.1 \
#		   src/external/libarchive-2.8.4])


# This adds the option of not compiling the pybeep library,
AC_ARG_ENABLE(pybeep,
  [  --disable-pybeep   Disable compilation of pybeep],
  [case "${enableval}" in
     yes | no ) WITH_PYBEEP="${enableval}" ;;
     *) AC_MSG_ERROR(bad value ${enableval} for --disable-pybeep) ;;
   esac],
  [WITH_PYBEEP="yes"]
)

dnl Make sure we register this option with Automake, so we know whether to
dnl descend into pybeep folder
AM_CONDITIONAL([WITH_PYBEEP], [test "x$WITH_PYBEEP" = "xyes"])

# This adds the option of not compiling the charm++ based parallel versions
AC_ARG_ENABLE(parallel,
  [  --disable-parallel   Disable compilation of parallel versions with charm++],
  [case "${enableval}" in
     yes | no ) WITH_CHARM_PLUS_PLUS="${enableval}" ;;
     *) AC_MSG_ERROR(bad value ${enableval} for --disable-parallel) ;;
   esac],
  [WITH_CHARM_PLUS_PLUS="yes"]
)

dnl Make sure we register this option with Automake, so we know whether to
dnl descend into charm++ folders
AM_CONDITIONAL([WITH_CHARM_PLUS_PLUS], [test "x$WITH_CHARM_PLUS_PLUS" = "xyes"])

# Option to forcibly disable OpenCL, even if it's available
AC_ARG_ENABLE(opencl,
  [  --disable-opencl   Don't use OpenCL even if it's available],
  [case "${enableval}" in
     yes | no ) IS_OK_TO_USE_OPENCL="${enableval}" ;;
     *) AC_MSG_ERROR(bad value ${enableval} for --disable-opencl) ;;
   esac],
  [IS_OK_TO_USE_OPENCL="yes"]
)
AM_CONDITIONAL([IS_OK_TO_USE_OPENCL], [test "x$IS_OK_TO_USE_OPENCL" = "xyes"])

# now we know it is (or is not) ok to use opencl, try and see if it is even available
AC_CONFIG_MACRO_DIR([m4])
AM_COND_IF([IS_OK_TO_USE_OPENCL], [AX_CHECK_CL])
AM_CONDITIONAL([WITH_OPENCL],[test "x$no_cl" = "x" -a "x$IS_OK_TO_USE_OPENCL" = "xyes"])

AC_OUTPUT
