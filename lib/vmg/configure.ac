
#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.60])
AC_INIT([VMG], [0.1], [tremolo-devel@ins.uni-bonn.de])

AC_MSG_NOTICE([****************************************************************])
AC_MSG_NOTICE([*      Configuring solver VMG                                  *])
AC_MSG_NOTICE([****************************************************************])

AC_CONFIG_SRCDIR([src/mg.cpp])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([1.7 -Wall foreign])

# We use mostly C++, but may need Fortran compiler infos for Lapack.
AC_LANG([C++])

# switch for one-sided communications
AC_MSG_CHECKING(whether to enable one-sided MPI2 communications)
AC_ARG_ENABLE([one_sided], 
	AS_HELP_STRING([--enable-one-sided], [Enable use of one-sided communications as according to MPI2 standard]),
	enable_one_sided=$enableval,
        enable_one_sided="no")
AC_MSG_RESULT($enable_one_sided)
AS_IF([test "$enable_one_sided" = "yes"], [
	    AC_DEFINE([VMG_ONE_SIDED], [1], [Enable use of one-sided MPI2 communications])
])

# several debug output switches
AC_MSG_CHECKING(whether to enable debugging)
AC_ARG_ENABLE([debug], 
	AS_HELP_STRING([--enable-debug], [Enable debugging]),
	enable_debug=$enableval,
	enable_debug="no")
AC_MSG_RESULT($enable_debug)
AS_IF([test "$enable_debug" = "yes"], [
	    CXXFLAGS=${CXXFLAGS-"-g -O0 -Wall"}
	    AC_DEFINE([DEBUG], [1], [Debugging output])
])

AC_MSG_CHECKING(whether to enable debug checking of matrices)
AC_ARG_ENABLE([debug-check-matrices], 
	AS_HELP_STRING([--enable-debug-check-matrices], [Enable matrix validity checks]),
	enable_debug_check_matrices=$enableval,
	enable_debug_check_matrices="no")
AC_MSG_RESULT($enable_debug_check_matrices)
AS_IF([test "$enable_debug_check_matrices" = "yes"], [
	    AC_DEFINE([DEBUG_MATRIX_CHECKS], [1], [Matrix validity checks])
])

AC_MSG_CHECKING(whether to enable debug output)
AC_ARG_ENABLE([debug-output], 
	AS_HELP_STRING([--enable-debug-output], [Enable debugging output]),
	enable_debug_output=$enableval,
	enable_debug_output="no")
AC_MSG_RESULT($enable_debug_output)
AS_IF([test "$enable_debug_output" = "yes"], [
	    AC_DEFINE([DEBUG_OUTPUT], [1], [Debugging output])
])

AC_MSG_CHECKING(whether to enable MPE profiling)
AC_ARG_ENABLE([mpe], 
	AS_HELP_STRING([--enable-mpe], [Enable MPE parallel profiling]),
	enable_mpe=$enableval,
	enable_mpe="no")
AC_MSG_RESULT($enable_mpe)

AC_MSG_CHECKING(whether to enable runtime measurement)
AC_ARG_ENABLE([timer], 
	AS_HELP_STRING([--enable-timer], [Enable runtime measurement]),
	enable_timer=$enableval,
	enable_timer="no")
AC_MSG_RESULT($enable_timer)
AS_IF([test "$enable_timer" = "yes"], [
	    AC_DEFINE([DEBUG_MEASURE_TIME], [1], [Runtime measurement])
])

AC_MSG_CHECKING(whether to enable runtime measurement output)
AC_ARG_ENABLE([timer-output],
        AS_HELP_STRING([--enable-timer-output], [Enable runtime measurement output]),
        enable_timer_output=$enableval,
        enable_timer_output="no")
AC_MSG_RESULT($enable_timer_output)
AS_IF([test "$enable_timer_output" = "yes"], [
            AC_DEFINE([DEBUG_MEASURE_TIME_OUTPUT], [1], [Output of runtime measurement])
])

AC_MSG_CHECKING(whether to enable MPI_Barrier debugging)
AC_ARG_ENABLE([debug-barrier], 
	AS_HELP_STRING([--enable-debug-barrier], [Enable MPI_Barrier debugging]),
	enable_debug_barrier=$enableval,
	enable_debug_barrier="no")
AC_MSG_RESULT($enable_debug_barrier)
AS_IF([test "$enable_debug_barrier" = "yes"], [
	    AC_DEFINE([DEBUG_BARRIER], [1], [MPI_Barrier debugging])
])

AC_MSG_CHECKING(whether to enable release version)
AC_ARG_ENABLE([release], 
	AS_HELP_STRING([--enable-release], [Compile with release settings]),
	enable_release=$enableval,
	enable_release="no")
AC_MSG_RESULT($enable_release)
AS_IF([test "$enable_release" = "yes"], [
	    AC_DEFINE([RELEASE], [1], [Release])
	    AC_DEFINE([NDEBUG], [1], [Skip asserts])
])

AC_ARG_VAR(BSPLINE_DEG, [Degree of interpolating B-Splines. Must be [3].])

if test -z "$BSPLINE_DEG"
then
	BSPLINE_DEG=3
fi

AC_MSG_CHECKING(interpolation b-spline degree)
AC_DEFINE_UNQUOTED([BSPLINE_DEGREE], ${BSPLINE_DEG}, [Interpolating B-Spline degree])
AC_MSG_RESULT(${BSPLINE_DEG})

# Checks for programs.
AM_MISSING_PROG([DOXYGEN], [doxygen])

# Get the C MPI compiler, because the following AC_FC_WRAPPERS macro works better if CC is from the same "family" as FC (which is MPI).
AX_PROG_CC_MPI(,,AC_MSG_FAILURE([The VMG solver requires an MPI C compiler.]))

# Do not even look for a non-MPI C++ compiler if MPICXX is set.
if test -n "$MPICXX" && test -z "$CXX"; then
  CXX=$MPICXX
fi
ACX_MPI([AC_DEFINE([HAVE_MPI], [1],
		   [Define if building the MPI parallel version])
	 CXX=$MPICXX
	 LIBS="$MPILIBS $LIBS"
	 AC_DEFINE([MPICH_SKIP_MPICXX], [1], [Skip C++ bindings])
	 AC_DEFINE([OMPI_SKIP_MPICXX], [1], [Skip C++ bindings])
	 AC_DEFINE([MPI_NO_CPPBIND], [1], [Skip C++ bindings])
	 AC_DEFINE([_MPICC_H], [1], [Skip C++ bindings])
	 AC_DEFINE([MPIBULL_SKIP_MPICXX], [1], [Skip C++ bindings])]
	 if test "$enable_mpe" = "yes"; then
	    AC_DEFINE([HAVE_MPE], [1], [Use MPE parallel profiler])
	 fi,
	[AC_DEFINE([sequent], [1],
		   [Define if building the serial version])])

# Get the Fortran MPI compiler, for Lapack.
AC_LANG_PUSH([Fortran])
# Do not even look for a non-MPI Fortran compiler if MPIFC is set.
AX_PROG_FC_MPI(,,AC_MSG_FAILURE([The VMG solver requires an MPI Fortran compiler.]))

# Find out how to call Fortran functions from C.
AC_FC_WRAPPERS
AC_LANG_POP([Fortran])

# Check for the restrict keyword.
AC_C_RESTRICT

AM_OPTIONS_VTK
AM_PATH_VTK([5.8.0],
	[AC_DEFINE([HAVE_VTK], [1], [VTK present on system])])

AC_SUBST([VTK_CXXFLAGS])
AC_SUBST([VTK_LDFLAGS])
AC_SUBST([VTK_LIBS])

AX_BOOST_BASE([1.34],[have_boost=yes],[have_boost=no])
if test "x$have_boost" = xyes ; then
  AX_BOOST_SYSTEM
  AX_BOOST_FILESYSTEM
  AX_BOOST_UNIT_TEST_FRAMEWORK
fi
AM_CONDITIONAL([COND_HAVE_BOOST_UNIT_TEST], [test "x$have_boost" = xyes && test x"$ax_cv_boost_unit_test_framework" = xyes])
AM_CONDITIONAL([COND_HAVE_BOOST_UNIT_TEST_LIB], [test "x$have_boost" = xyes && test x"$ax_cv_boost_unit_test_framework" = xyes && test x"$ax_cv_boost_unit_test_framework_lib" = xyes])

AC_PROG_RANLIB

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])
AX_LAPACK
AM_CONDITIONAL([HAVE_LAPACK], [test "x$ax_lapack_ok" = xyes])

# Checks for header files.

# Checks for types.

# Set up FCS library arguments.
m4_ifdef([AX_FCS_CONFIG_SOLVER_ARGS],[AX_FCS_CONFIG_SOLVER_ARGS])

# Set up FCS library.
m4_ifdef([AX_FCS_CONFIG_SOLVER],[AX_FCS_CONFIG_SOLVER])

# Set up FCS types 
m4_ifdef([AX_FCS_TYPES],
  [AC_DEFINE([vmg_float], [fcs_float],
     [Define to the primary interface and computation floating type.])
   AC_DEFINE([vmg_int], [fcs_int],
     [Define to the primary integer type for particle indices.])],
  [AC_DEFINE([vmg_float], [double],
     [Define to the primary interface and computation floating type.])
   AC_DEFINE([vmg_int], [int],
     [Define to the primary integer type for particle indices.])])

# Set up information for FCS package if inside the FCS tree.
m4_ifdef([AX_FCS_PACKAGE_RESET],[
AX_FCS_PACKAGE_RESET
AX_FCS_PACKAGE_ADD([vmg_LIBS],[-lfcs_vmg])
AX_FCS_PACKAGE_ADD([vmg_LIBS_A],[src/libfcs_vmg.a])
AX_FCS_PACKAGE_ADD([vmg_LDADD],[$VTK_LDFLAGS $VTK_LIBS $BOOST_SYSTEM_LIB $BOOST_FILESYSTEM_LIB $LAPACK_LIBS $BLAS_LIBS])
AX_FCS_PACKAGE_ADD([CXXLIBS_USE],[yes])
])

# Checks for structures.

# Checks for compiler characteristics.

# Checks for library functions.

# Checks for system services

# we save the cache here because the tests below fail more commonly

AC_CONFIG_FILES([Makefile
	src/Makefile
	doc/Makefile
	test/Makefile])
AC_OUTPUT
