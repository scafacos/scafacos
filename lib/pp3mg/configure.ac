# Initialize Autoconf.  We require version 2.59 or newer.
AC_PREREQ([2.59])
AC_INIT([pmg], [1.0])

AC_MSG_NOTICE([****************************************************************])
AC_MSG_NOTICE([*      Configuring solver PP3MG                                *])
AC_MSG_NOTICE([****************************************************************])

# Share macro files with the rest of scafacos.
AC_CONFIG_MACRO_DIR([m4])
# Share helper scripts with the rest of scafacos.
AC_CONFIG_AUX_DIR([build-aux])

# Initialize Automake; require 1.9.6 or newer,
# do not enforce GNU rules.
AM_INIT_AUTOMAKE([1.9.6 foreign])

AX_PROG_CC_MPI(,,AC_MSG_FAILURE([The PP3MG solver requires an MPI C compiler.]))

# The AC_PROG_CC_C99 macro is relatively new in Autoconf, so
# if the developer only used an older version, the user will
# have to specify the needed flags (e.g., CC="gcc -std=gnu99").
m4_ifdef([AC_PROG_CC_C99], [AC_PROG_CC_C99])
# Do we -lm?
AC_SEARCH_LIBS([ceil], [m])
# We need the C99 round function.
AC_CHECK_DECL([round], [],
  [AC_MSG_ERROR([need declaration of C99 round function])],
  [AC_INCLUDES_DEFAULT
   [#include <math.h>]])

# For Fortran, we also just use the MPIFC throughout.
AC_LANG_PUSH([Fortran])
# Do not even look for a non-MPI Fortran compiler if MPIFC is set.
AX_PROG_FC_MPI(,,AC_MSG_FAILURE([The PP3MG solver requires an MPI Fortran compiler.]))

AC_FC_FREEFORM

# Accept preprocessed source files names ending in ".f90".
AC_FC_PP_SRCEXT([f90])

# Get the flag to use for setting Fortran defines on the command line.
AC_FC_PP_DEFINE

# Set FCDEFS and adjust -D flag if needed.
FCDEFS="$DEFS -DHAVE_FCONFIG_H"
if test "x$FC_DEFINE" != x-D; then
  FCDEFS=`echo " $FCDEFS" | sed "s@ -D@ $FC_DEFINE@g"`
fi
AC_SUBST([FCDEFS])

# It is probably confusing for the user to have both CPPFLAGS and FCCPPFLAGS
# (but it may be necessary to set them differently in some situations).  Since
# we don't use any preprocessor tests below, default the latter to the former
# at 'make' run time if they are equal, so setting the former affects both.
if test "x$CPPFLAGS" = "x$FCCPPFLAGS" || test -z "$FCCPPFLAGS"; then
  FCCPPFLAGS=`echo " $CPPFLAGS" | sed "s@ -D@ $FC_DEFINE@g"`
fi
AC_SUBST([FCCPPFLAGS])

# Find out how to add to the Fortran include path.
AC_FC_MODULE_FLAG
# Find out the Fortran module file extension.
AC_FC_MODULE_EXTENSION
AC_LANG_POP([Fortran])

# Checks for typedefs, structures, and compiler characteristics.

# Init libtool
LT_INIT([disable-shared])

######################################################################
# Set up parameters.
######################################################################

# Set up FCS library arguments.
m4_ifdef([AX_FCS_CONFIG_SOLVER_ARGS],[AX_FCS_CONFIG_SOLVER_ARGS])

######################################################################
# Set up FCS and solver stuff.
######################################################################

# Set up FCS library.
m4_ifdef([AX_FCS_CONFIG_SOLVER],[AX_FCS_CONFIG_SOLVER])

# Set up information for FCS package if inside the FCS tree.
m4_ifdef([AX_FCS_PACKAGE_RESET],[
AX_FCS_PACKAGE_RESET
AX_FCS_PACKAGE_ADD([pp3mg_LIBS],[-lfcs_pp3mg])
AX_FCS_PACKAGE_ADD([pp3mg_LIBS_A],[libfcs_pp3mg.la])
])

# Create output files and headers.
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile pmg/Makefile pp3mg/Makefile])
AC_OUTPUT
