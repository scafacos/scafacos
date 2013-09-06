# Initialize Autoconf.  We require version 2.59 or newer.
AC_PREREQ([2.59])
AC_INIT([PEPC], [rev3965], [pepc@fz-juelich.de], [pepc])

AC_MSG_NOTICE([****************************************************************])
AC_MSG_NOTICE([*      Configuring solver PEPC                                 *])
AC_MSG_NOTICE([****************************************************************])

# Share macro files with the rest of scafacos.
AC_CONFIG_MACRO_DIR([m4])
# Share helper scripts with the rest of scafacos.
AC_CONFIG_AUX_DIR([build-aux])

# Initialize Automake; require 1.9.6 or newer, do not enforce GNU rules.
AM_INIT_AUTOMAKE([1.9.6 foreign])

# Checks for programs.

# Get the C MPI compiler.
AX_PROG_CC_MPI(,,AC_MSG_FAILURE([The PEPC solver requires an MPI C compiler.]))
   
# Find out whether the C compiler accepts '-c -o'. (What's this good for?)
AM_PROG_CC_C_O

# Get the Preprocessor.
AC_PROG_CPP

# Get the Fortran MPI compiler.
AX_PROG_FC_MPI(,,AC_MSG_FAILURE([The PEPC solver requires an MPI Fortran compiler.]))

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

# Get the Fortran compiler to accept free-form source.
AC_FC_FREEFORM
# Get it to accept long source file lines.
# Unset the cache variable to force a detection with the currently selected fixed- or free-form.
unset ac_cv_fc_line_length
AC_FC_LINE_LENGTH([unlimited])
# Find out how to add to the Fortran include path.
AC_FC_MODULE_FLAG
# Find out the Fortran module file extension.
AC_FC_MODULE_EXTENSION

AC_PROG_RANLIB

AC_CHECK_TOOLS([AR], [ar xiar], [false])
: ${ARFLAGS=r}
AC_SUBST([ARFLAGS])

: ${RM="rm -f"}
AC_SUBST([RM])

# for some code parts, this solver needs OpenMP, so check whether its available
# this macro sets OPENMP_CFLAGS, OPENMP_CXXFLAGS, OPENMP_FCFLAGS, depending on choice of current langugage - we need the Fortran here
AC_LANG_PUSH([Fortran])dnl
AX_OPENMP
AC_LANG_POP([Fortran])dnl

######################################################################
# Set up parameters.
######################################################################

# Set up FCS library arguments.
AX_FCS_CONFIG_SOLVER_ARGS

######################################################################
# Checks for libraries.
######################################################################

# Set up the sorting library.
AX_SORT_LIB

######################################################################
# Set up FCS and solver stuff.
######################################################################

# Set up FCS library.
AX_FCS_CONFIG_SOLVER

# we need OpenMP for some parts of the solver
FCFLAGS="$FCFLAGS $OPENMP_FCFLAGS"

# configure sub-library
AC_CONFIG_SUBDIRS([openpa])

# Set up information for FCS package if inside the FCS tree.
m4_ifdef([AX_FCS_PACKAGE_RESET],[
AX_FCS_PACKAGE_RESET
AX_FCS_PACKAGE_ADD([pepc_LIBS],[-lfcs_pepc])
AX_FCS_PACKAGE_ADD([pepc_LIBS_A],[libfcs_pepc.a])
AX_FCS_PACKAGE_ADD([FCLIBS_USE],[yes])
AX_FCS_PACKAGE_ADD([FCOMP_USE],[yes])
])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_SIZEOF([int *])

# Find out the default type sizes.
AX_FORTRAN_KINDS_DEFAULT

# Checks for library functions.

# Create output files.

# Remove obsolete module dependency rules (for convenience only, REMOVE later).
m4_esyscmd([rm -f src/module.deps])

# These are from automake-generated .in files.
AC_CONFIG_FILES([Makefile
		 src/Makefile
		 sl_pepc/Makefile])

# Create config headers: config.h comes from autoheader-generated config.h.in,
# fconfig.h.in is written manually and only contains stuff needed for Fortran.
AC_CONFIG_HEADERS([config.h])
# Convert config.h into an fconfig.h suitable for inclusion into Fortran.
# The sed script transforms C style /* comments */ into Fortran style.
AC_CONFIG_HEADERS([fconfig.h],
 [[sed '
    /\/\*\(.*\)\*\//{
      s,,!\1,p
      d
    }
    /\/\*/,/\*\//{
      s,/\*,,
      s,\*/,,
      s,^,!,p
      d
    }' fconfig.h > fconfig.ht
   mv -f fconfig.ht fconfig.h
]])

AC_OUTPUT
