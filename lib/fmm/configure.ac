# Initialize Autoconf.  We require version 2.59 or newer.
AC_PREREQ([2.59])
AC_INIT([FMM], [0.1], [libfm@fz-juelich.de])

AC_MSG_NOTICE([****************************************************************])
AC_MSG_NOTICE([*      Configuring solver FMM                                  *])
AC_MSG_NOTICE([****************************************************************])

# Share macro files with the rest of scafacos.
AC_CONFIG_MACRO_DIR([m4])
# Share helper scripts with the rest of scafacos.
AC_CONFIG_AUX_DIR([build-aux])

# Initialize Automake; require 1.11 or newer, do not enforce GNU rules.
# We requires version 1.11 so subdir-objects work correctly with Fortran.
AM_INIT_AUTOMAKE([1.11 foreign subdir-objects])

######################################################################
# Checks for programs. 
######################################################################

# Get the C MPI compiler.
AX_PROG_CC_MPI(,,AC_MSG_FAILURE([The FMM solver requires an MPI C compiler.]))

# Find out whether the C compiler accepts '-c -o'. (What's this good for?)
AM_PROG_CC_C_O

# Get the Fortran MPI compiler.
AX_PROG_FC_MPI(,,AC_MSG_FAILURE([The FMM solver requires an MPI Fortran compiler.]))
AC_SUBST([F77], ['$(FC)'])

# Get the Preprocessor. FC variable should be set
ACX_PROG_FCCPP

# Check whether gfortran >= 4.6 is used and add option --fno-whole-file
$FC --version >conftest.log 2>&1
gfortran_id_string=`grep 'GNU Fortran' conftest.log | head -1`
if test -n "$gfortran_id_string" ; then
  [gfortran_version=`echo "$gfortran_id_string" | sed -e 's/.*\([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/'`]
  [gfortran_version2=`echo "$gfortran_version" | sed -e 's/\.//g'`]
  while test "${#gfortran_version2}" -lt "3" ; do
    gfortran_version2="${gfortran_version2}0"
  done
  if test "$gfortran_version2" -ge "460" ; then
    case $FCFLAGS in
      *-fno-whole-file*) ;;
      *)
        AC_MSG_NOTICE([found gfortran version $gfortran_version, enabling flag -fno-whole-file])
        FCFLAGS="$FCFLAGS -fno-whole-file"
        ;;
    esac
  fi
fi

# Accept preprocessed source files names ending in ".f90" and ".f".
AC_FC_PP_SRCEXT([f90])
AC_FC_PP_SRCEXT([f])

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
save_FCFLAGS=$FCFLAGS
AC_FC_FREEFORM
FCFLAGS=$save_FCFLAGS
case $ac_cv_fc_freeform in
  none | unknown) FC_FREEFORM= ;;
  *) FC_FREEFORM=$ac_cv_fc_freeform ;;
esac
AC_SUBST([FC_FREEFORM])

# Get the Fortran compiler to accept fixed-form source, even when free-form
# is the default.
save_FCFLAGS=$FCFLAGS
FCFLAGS="$FCFLAGS $FC_FREEFORM"
AC_FC_FIXEDFORM
FCFLAGS=$save_FCFLAGS
case $ac_cv_fc_fixedform in
  none | unknown) FC_FIXEDFORM= ;;
  *) FC_FIXEDFORM=$ac_cv_fc_fixedform ;;
esac
AC_SUBST([FC_FIXEDFORM])

# Get it to accept long source file lines.
# Unset the cache variable to force a detection with the currently selected fixed- or free-form.
unset ac_cv_fc_line_length
AC_FC_LINE_LENGTH

# Find out how to turn off range checking, but don't turn it off by default.
save_FCFLAGS=$FCFLAGS
AX_FC_NO_RANGE_CHECK
FCFLAGS=$save_FCFLAGS
case $ax_cv_fc_no_range_check in
  none | unknown) FC_NO_RANGE_CHECK= ;;
  *) FC_NO_RANGE_CHECK=$ax_cv_fc_no_range_check ;;
esac
AC_SUBST([FC_NO_RANGE_CHECK])

# Find out how to set the default integer kind to 8 bytes.
AX_FC_INTEGER_KIND([8])
case $ax_cv_fc_integer_kind_8 in
  none | unknown) FC_INTEGER_SIZE_8= ;;
  *) FC_INTEGER_SIZE_8=$ax_cv_fc_integer_kind_8 ;;
esac
AC_SUBST([FC_INTEGER_SIZE_8])

# Find out how to add to the Fortran include path.
AC_FC_MODULE_FLAG
# Find out the file name extensions of module files.
AC_FC_MODULE_EXTENSION

AC_CHECK_TOOLS([AR], [ar xiar], [false])
: ${ARFLAGS=r}
AC_SUBST([ARFLAGS])

: ${RM="rm -f"}
AC_SUBST([RM])

# Init libtool
LT_INIT([disable-shared])

######################################################################
# Set up parameters.
######################################################################

# Set up FCS library arguments.
AX_FCS_CONFIG_SOLVER_ARGS

# Set up FMM solver arguments.
AX_FCS_FMM_ARGS

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

# Set up FMM solver.
AX_FCS_FMM_SOLVER

# Find out whether intrinsics are available.
AX_CHECK_IBM_INTRINSICS
AX_CHECK_SSE2_INTRINSICS

AX_COMPILER_VENDOR
AC_CHECK_SIZEOF([fcs_float])

if test "x${ax_intrinsics_ibm_f}" = xyes -a "x${ac_cv_sizeof_fcs_float}" = x8 ; then
  fmm_ppc450_64bit=yes
fi
if test "x${ax_intrinsics_sse2_c}" = xyes -a "x${ac_cv_sizeof_fcs_float}" = x8 ; then
  if test "x${ax_cv_c_compiler_vendor}" = xintel ; then
    fmm_sse2_64bit=yes
  fi
fi

AM_CONDITIONAL([ENABLE_IBM_F_INTRINSICS],[test "x${ax_intrinsics_ibm_f}" = xyes])

#AM_CONDITIONAL([ENABLE_PPC450_32BIT],[test "x${fmm_ppc450_32bit}" = xyes])
AM_CONDITIONAL([ENABLE_PPC450_64BIT],[test "x${fmm_ppc450_64bit}" = xyes])
#AM_CONDITIONAL([ENABLE_SSE2_32BIT],[test "x${fmm_sse2_32bit}" = xyes])
AM_CONDITIONAL([ENABLE_SSE2_64BIT],[test "x${fmm_sse2_64bit}" = xyes])
AM_CONDITIONAL([ENABLE_GENERIC],[test "x${fmm_ppc450_32bit}${fmm_ppc450_64bit}${fmm_sse2_32bit}${fmm_sse2_64bit}" != xyes])

# Add FMM solver stuff to main FCS package.
AX_FCS_PACKAGE_RESET
AX_FCS_PACKAGE_ADD([fmm_LIBS],[-lfcs_fmm])
AX_FCS_PACKAGE_ADD([fmm_LIBS_A],[libfcs_fmm.la])
if test "x${FMM_MP}" = xFMM_MP_ARMCI ; then
  AX_FCS_PACKAGE_ADD([SUB_PACKAGES],[armci])
fi
if test "x${FMM_MP}" = xFMM_MP_A1 ; then
  AX_FCS_PACKAGE_ADD([SUB_PACKAGES],[a1])
fi
if test "x${FMM_MP}" = xFMM_MP_ARMCIMPI ; then
  AX_FCS_PACKAGE_ADD([SUB_PACKAGES],[armci-mpi])
fi
# The following is currently not required, because -lpami is automatically detected as one of the required MPI libraries.
#if test "x${FMM_MP}" = xFMM_MP_SIMPLE_ARMCI -a "x${enable_simple_armci_device}" = xpamid ; then
#  AX_FCS_PACKAGE_ADD([fmm_LIBS],[-lpami])
#fi
AX_FCS_PACKAGE_ADD([FCLIBS_USE],[yes])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

######################################################################
# Create output files.
######################################################################

# These are from automake-generated .in files.
AC_CONFIG_FILES([Makefile
  src/Makefile
  sl_fmm/Makefile])

#if test "x${FMM_MP}" = xFMM_MP_SIMPLE_ARMCI ; then
#  AC_CONFIG_FILES([simple-armci/Makefile
#    simple-armci/generic/Makefile
#    simple-armci/dcmfd/Makefile])
#fi

m4_foreach_w([mm], [01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50],
 [AM_CONDITIONAL([ENABLE_FMM_MM_]mm, [test "mm" -le "$enable_fcs_fmm_max_mpol"])
 ])

AX_FILELIST([filelist.am],[src/unrolled/fmmgradt],[*.f],[*main.f *p00.f],,[ENABLE_FMM_MM_`echo ${ax_filelist_name} | sed -e 's/.*p//' -e 's/\..*//'`])
AC_CONFIG_FILES([src/unrolled/fmmgradt/Makefile])

AX_FILELIST([filelist.am],[src/unrolled/fmmmopn],[*.f],[*main.f *p00.f],,[ENABLE_FMM_MM_`echo ${ax_filelist_name} | sed -e 's/.*p//' -e 's/\..*//'`])
AC_CONFIG_FILES([src/unrolled/fmmmopn/Makefile])

AX_FILELIST([filelist.am],[src/unrolled/fmmoopn],[*.f],[*main.f *p00.f],,[ENABLE_FMM_MM_`echo ${ax_filelist_name} | sed -e 's/.*p//' -e 's/\..*//'`])
AC_CONFIG_FILES([src/unrolled/fmmoopn/Makefile])

# Generate file lists for "errorcontrol" directories.
m4_foreach_w([dir], [errorcontrol/0d errorcontrol/1d errorcontrol/2d errorcontrol/3d],
 [AX_FILELIST([filelist.am],[src/]dir,[*.h *.c *.f *.f90 include/*.h])[]dnl
  AC_CONFIG_FILES([src/]dir[/Makefile])[]dnl
 ])

# Remove obsolete module dependency rules (for convenience only, REMOVE later).
m4_esyscmd([rm -f src/module.deps])

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
