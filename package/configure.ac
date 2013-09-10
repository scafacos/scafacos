AC_PREREQ([2.59])
AC_INIT([ScaFaCoS FCS], [0.1], [scafacos-dvlp@fz-juelich.de], [scafacos-fcs],
  [http://www2.fz-juelich.de/jsc/scafacos/])
m4_ifndef([AC_PACKAGE_URL],
	  [AC_SUBST([PACKAGE_URL], [http://www.fz-juelich.de/jsc/scafacos/])])

AC_MSG_NOTICE([****************************************************************])
AC_MSG_NOTICE([*      Configuring ScaFaCos package                            *])
AC_MSG_NOTICE([****************************************************************])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([build-aux])

AM_INIT_AUTOMAKE([1.9.6 foreign])

# Print an error if a (custom) AX_ or ACX_ macro appears in the ouput
m4_pattern_forbid([^A[C]?X_])

# Get the C MPI compiler.
AX_PROG_CC_MPI(,,AC_MSG_FAILURE([ScaFaCoS requires an MPI C compiler.]))
# Get the C++ MPI compiler.
AX_PROG_CXX_MPI(,,AC_MSG_FAILURE([ScaFaCoS requires an MPI C++ compiler.]))
# Get the Fortran MPI compiler.
AX_PROG_FC_MPI(,,AC_MSG_FAILURE([ScaFaCoS requires an MPI Fortran compiler.]))

AC_PROG_RANLIB

# Add single-lib feature
AX_FCS_ENABLE_SINGLE_LIB_ARG
AX_FCS_ENABLE_SINGLE_LIB

# List of all solvers.
m4_define([all_solver_methods], [direct ewald fmm memd mmm1d mmm2d p2nfft p3m pepc pp3mg vmg wolf])

m4_define([all_common_methods], [near gridsort resort fcs_common pnfft pfft fftw3_mpi_common fftw3_common fftw3_mpi fftw3])

# Determine enabled solvers.
AC_ARG_ENABLE([fcs-package-solvers],[],[],[enable_fcs_package_solvers="all_solver_methods"])

#echo "enable_fcs_package_solvers: $enable_fcs_package_solvers"

save_IFS=$IFS
IFS=", "
for s in ${enable_fcs_package_solvers} ; do
  enabled_solvers="${enabled_solvers} ${s} "
done
IFS=$save_IFS

#echo "enabled_solvers: |$enabled_solvers|"

function process_package
{
  local topdir="$1"
  local pkgdir="$2"

#  echo "processing: $pkgdir"

  AX_FCS_PACKAGE_PROCESS([${topdir}],[$pkgdir],[fcs ]all_solver_methods all_common_methods)

#  AC_FOREACH([method],[fcs ]all_solver_methods all_common_methods,[
#    echo "method:"
#    echo "  USE:       $ax_fcs_package_[]method[]_USE"
#    echo "  LIBS:      $ax_fcs_package_[]method[]_LIBS"
#    echo "  LIBS_A:    $ax_fcs_package_[]method[]_LIBS_A"
#    echo "  LDADD:     $ax_fcs_package_[]method[]_LDADD"
#  ])
}

# process main package information
process_package ".." "."

# process solver package information
AC_FOREACH([solver], all_solver_methods,[
  if test "x${enabled_solvers%* solver *}" != "x${enabled_solvers}" ; then
    process_package ".." "lib/solver"
    ax_fcs_package_[]solver[]_USE=yes
  fi
])

if test "x$enable_single_lib" = xyes ; then
  SCAFACOS_PC_LIBS="${SCAFACOS_PC_LIBS} -lfcs"
  SCAFACOS_MK_LDADD="${SCAFACOS_MK_LDADD} \$(top_builddir)/package/libfcs.a"
  SCAFACOS_MK_DEPS="${SCAFACOS_MK_DEPS} \$(top_builddir)/package/libfcs.a"
fi

if test "x${ax_fcs_package_COMP_USE}" != x ; then
  AC_LANG_PUSH([C])
  AX_OPENMP
  AC_LANG_POP([C])
  CFLAGS="$CFLAGS $OPENMP_CFLAGS"
fi
if test "x${ax_fcs_package_FCOMP_USE}" != x ; then
  AC_LANG_PUSH([Fortran])
  AX_OPENMP
  AC_LANG_POP([Fortran])
  FCFLAGS="$FCFLAGS $OPENMP_FCFLAGS"
fi
if test "x${ax_fcs_package_CXXOMP_USE}" != x ; then
  AC_LANG_PUSH([C++])
  AX_OPENMP
  AC_LANG_POP([C++])
  CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
fi

AC_FOREACH([mod],[fcs ]all_solver_methods all_common_methods,[
if test "x${ax_fcs_package_[]mod[]_USE}" != x ; then
  if test "x$enable_single_lib" != xyes ; then
    SCAFACOS_PC_LIBS="${SCAFACOS_PC_LIBS} ${ax_fcs_package_[]mod[]_LIBS}"
  fi
  SCAFACOS_PC_LIBS="${SCAFACOS_PC_LIBS} ${ax_fcs_package_[]mod[]_LDADD}"
  LIBS_A=
  LIBS_A_PKG=
  for d in ${ax_fcs_package_[]mod[]_LIBS_A} ; do
    if test "x${d#/}" = "x${d}" ; then
      LIBS_A="${LIBS_A} \$(top_builddir)/${d}"
      LIBS_A_PKG="${LIBS_A_PKG} \$(top_builddir)/../${d}"
    else
      LIBS_A="${LIBS_A} ${d}"
      LIBS_A_PKG="${LIBS_A_PKG} ${d}"
    fi
  done
  if test "x$enable_single_lib" != xyes ; then
    SCAFACOS_MK_LDADD="${SCAFACOS_MK_LDADD} ${LIBS_A}"
    SCAFACOS_MK_DEPS="${SCAFACOS_MK_DEPS} ${LIBS_A}"
  fi
  SCAFACOS_MK_LDADD="${SCAFACOS_MK_LDADD} ${ax_fcs_package_[]mod[]_LDADD}"
  SCAFACOS_MK_LIBS_A="${SCAFACOS_MK_LIBS_A} ${LIBS_A_PKG}"
fi
])

# It's ugly, but probably safer: add full Fortran or C++ libraries if required.
use_fcs_full_ldflags=yes

# Unset the cache variables for Fortran and C++ libraries (before the following checks),
# because the compiler flags here may be different from the top-level configure (e.g., OpenMP flags!).

# Determine which Fortran libraries are required by a non-Fortran linker.
if test "x${ax_fcs_package_FCLIBS_USE}" != x ; then
  AS_UNSET([ac_cv_fc_libs])
  if test "x${use_fcs_full_ldflags}" = xyes ; then
    AC_FC_LIBRARY_LDFLAGS
    m4_ifdef([ACX_FC_LIBRARY_LDFLAGS_FIX],[ACX_FC_LIBRARY_LDFLAGS_FIX])
    SCAFACOS_PC_LIBS="${SCAFACOS_PC_LIBS} ${FCLIBS}"
    SCAFACOS_MK_LDADD="${SCAFACOS_MK_LDADD} ${FCLIBS}"
  else
    AS_UNSET([ax_cv_diff_library_ldflags_fc])
    AX_DIFF_LIBRARY_LDFLAGS_FC
    SCAFACOS_PC_LIBS="${SCAFACOS_PC_LIBS} ${FCLIBS_DIFF}"
    SCAFACOS_MK_LDADD="${SCAFACOS_MK_LDADD} ${FCLIBS_DIFF}"
  fi
fi

# Determine which C++ libraries are required by a non-C++ linker.
if test "x${ax_fcs_package_CXXLIBS_USE}" != x ; then
  AS_UNSET([ax_cv_cxxlibs])
  if test "x${use_fcs_full_ldflags}" = xyes ; then
    AX_CXX_LIBRARY_LDFLAGS
    SCAFACOS_PC_LIBS="${SCAFACOS_PC_LIBS} ${CXXLIBS}"
    SCAFACOS_MK_LDADD="${SCAFACOS_MK_LDADD} ${CXXLIBS}"
  else
    AS_UNSET([ax_cv_diff_library_ldflags_cxx])
    AX_DIFF_LIBRARY_LDFLAGS_CXX
    SCAFACOS_PC_LIBS="${SCAFACOS_PC_LIBS} ${CXXLIBS_DIFF}"
    SCAFACOS_MK_LDADD="${SCAFACOS_MK_LDADD} ${CXXLIBS_DIFF}"
  fi
fi

# Remove leading, trailing and multiple white spaces.
z= ; for x in ${SCAFACOS_PC_LIBS} ; do z="${z} ${x}" ; done ; SCAFACOS_PC_LIBS="${z# }"
z= ; for x in ${SCAFACOS_MK_LDADD} ; do z="${z} ${x}" ; done ; SCAFACOS_MK_LDADD="${z# }"
z= ; for x in ${SCAFACOS_MK_DEPS} ; do z="${z} ${x}" ; done ; SCAFACOS_MK_DEPS="${z# }"
z= ; for x in ${SCAFACOS_MK_LIBS_A} ; do z="${z} ${x}" ; done ; SCAFACOS_MK_LIBS_A="${z# }"

#echo "SCAFACOS_PC_LIBS:   $SCAFACOS_PC_LIBS"
#echo "SCAFACOS_MK_LDADD:  $SCAFACOS_MK_LDADD"
#echo "SCAFACOS_MK_DEPS:   $SCAFACOS_MK_DEPS"
#echo "SCAFACOS_MK_LIBS_A: $SCAFACOS_MK_LIBS_A"

AC_SUBST([SCAFACOS_PC_LIBS])
AC_SUBST([SCAFACOS_MK_LDADD])
AC_SUBST([SCAFACOS_MK_DEPS])
AC_SUBST([SCAFACOS_MK_LIBS_A])

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([scafacos-fcs.pc])
AC_CONFIG_FILES([scafacos-fcs.mk])

AC_OUTPUT

AC_MSG_NOTICE([])
AC_MSG_NOTICE([****************************************************************])
AC_MSG_NOTICE([*                   Package summary                            *])
AC_MSG_NOTICE([****************************************************************])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([PC_LIBS  = ${SCAFACOS_PC_LIBS}])
AC_MSG_NOTICE([MK_LDADD = ${SCAFACOS_MK_LDADD}])
AC_MSG_NOTICE([MK_DEPS  = ${SCAFACOS_MK_DEPS}])
