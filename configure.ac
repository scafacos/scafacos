# Copyright (C) 2011 The ScaFaCoS project
#  
# This file is part of ScaFaCoS.
#  
# ScaFaCoS is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
#  ScaFaCoS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser Public License for more details.
#  
#  You should have received a copy of the GNU Lesser Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
#
# Initialize Autoconf.  We require version 2.59 or newer.
AC_PREREQ([2.59])
define([version_suffix],[])
AC_INIT([ScaFaCoS fcs], [0.1version_suffix], [scafacos-dvlp@fz-juelich.de], [scafacos-fcs],
	[http://www.fz-juelich.de/jsc/scafacos/])
m4_ifndef([AC_PACKAGE_URL],
	  [AC_SUBST([PACKAGE_URL], [http://www.fz-juelich.de/jsc/scafacos/])])

AC_MSG_NOTICE([****************************************************************])
AC_MSG_NOTICE([*      Configuring ScaFaCos main project                       *])
AC_MSG_NOTICE([****************************************************************])

# Additional macros are stored here.
AC_CONFIG_MACRO_DIR([m4])
# Additional helper scripts are stored here.
AC_CONFIG_AUX_DIR([build-aux])

# Initialize Automake; require 1.9.6 or newer, do not enforce GNU rules, create pax/posix tarfiles in make dist to support longer filenames in doxgen stuff.
AM_INIT_AUTOMAKE([1.9.6 -Wall -Wno-override foreign tar-pax])

# Print an error if a (custom) AX_ or ACX_ macro appears in the ouput
m4_pattern_forbid([^A[C]?X_])

# List of all solvers.
m4_define([all_solver_methods], [direct ewald fmm memd mmm1d mmm2d p2nfft p3m pepc pp3mg vmg wolf])

# List of solvers that have their own configure
m4_define([subconfigure_solver_methods], [fmm p2nfft pepc pp3mg vmg])

# Solvers which are licensed under GPL or depend on GPL software.
m4_define([gpl_solver_methods], [p2nfft p3m mmm1d mmm2d])

# Checks for programs.
# Get the C MPI compiler.
AX_PROG_CC_MPI(,,AC_MSG_FAILURE([ScaFaCoS requires an MPI C compiler.]))

# Get the Preprocessor.
AC_PROG_CPP

# Get the C++ MPI compiler.
AX_PROG_CXX_MPI(,,AC_MSG_FAILURE([ScaFaCoS requires an MPI C++ compiler.]))

# Get the Fortran MPI compiler.
AX_PROG_FC_MPI(,,AC_MSG_FAILURE([ScaFaCoS requires an MPI Fortran compiler.]))

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
AC_FC_LINE_LENGTH
# Find out how to add to the Fortran include path.
AC_FC_MODULE_FLAG
# Find out the Fortran module file extension.
AC_FC_MODULE_EXTENSION
# Find out how to link against Fortran libraries.
AC_FC_LIBRARY_LDFLAGS
# Bugfix for BlueGene xlf2003: remove spurious -link entry.
ACX_FC_LIBRARY_LDFLAGS_FIX
# Find out how to name a MAIN Fortran entry point in C.
AC_FC_MAIN
# Find out how to call Fortran functions from C.
AC_FC_WRAPPERS

# Find out whether we need "use ifport" for rand.
AC_CACHE_CHECK([whether to "use ifport" for rand],
 [fcs_cv_fc_module_ifport],
 [fcs_cv_fc_module_ifport=no
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
   [dnl mind Fortran indenting
      use ifport
      x = rand()])],
   [fcs_cv_fc_module_ifport=yes])
  AC_LANG_POP([Fortran])])
if test "x$fcs_cv_fc_module_ifport" = xyes ; then
  AC_DEFINE([NEED_USE_IFPORT], [1],
   [Define if 'use ifport' is needed for functions like 'rand'.])
fi

AC_PROG_RANLIB
AM_PROG_CC_C_O
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

# Checks for libraries.
AC_SEARCH_LIBS([log], [m])
AC_SEARCH_LIBS([gzopen], [z], [have_zlib=yes])

# Checks for header files.
AC_CHECK_HEADERS([zlib.h])

# check for iomanip: iomanip doesn't work in certain combinations of
# icc and gcc.
AC_LANG_PUSH(C++)
AC_CHECK_HEADERS([iomanip], [have_iomanip=yes], [have_iomanip=no])
AM_CONDITIONAL([HAVE_IOMANIP], [test x$have_iomanip = xyes])
AC_LANG_POP

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_SIZEOF([int *])
AC_PROG_CC_C99

# Checks for library functions.

#######################################################################
# FCS library section
#######################################################################

# Set up FCS library stuff.
AX_FCS_CONFIG_TOP

# Let's see if we can get 'make check' to start an MPI job automatically.
ACX_MPIEXEC

#######################################################################
# Solver section
#######################################################################

# If a only a single libfcs.a should be installed, then prevent install
# of sub-package libraries using option --disable-library-install.
if test "x$enable_single_lib" = xyes ; then
  case $ac_configure_args in
    *--disable-library-install*) ;;
    *) ac_configure_args="$ac_configure_args --disable-library-install" ;;
  esac
fi

# Decide which solvers we want to build.
AC_ARG_ENABLE([fcs-solvers],
 [AS_HELP_STRING([--enable-fcs-solvers=all|SOLVER,...],
  [choose solver methods to build @<:@all@:>@ (]all_solver_methods[)])],,
  [enable_fcs_solvers=all])
if test "x$enable_fcs_solvers" = xall ; then
  enable_fcs_solvers="all_solver_methods"
fi

#echo "enabled solvers: ${enable_fcs_solvers}"

# check all enabled solvers and check whether their source is available
save_IFS=$IFS
IFS=", "
for solver in $enable_fcs_solvers ; do
  case " all_solver_methods " in
    *\ $solver\ *)
      if test -d "$srcdir/lib/$solver"; then
        eval "use_fcs_$solver=yes"
      else
        AC_MSG_WARN([disabling solver $solver as it is not present in source tree])
        eval "use_fcs_$solver=\"$solver not present in source tree\""
      fi
      ;;
    *)
      AC_MSG_WARN([unknown solver: $solver])
      ;;
   esac
done
IFS=$save_IFS

# add option to disable GPL only solvers
AC_ARG_ENABLE([gpl],
 [AS_HELP_STRING([--disable-gpl],[disable GPL solvers])],,
 [enable_gpl=yes])
if test "x$enable_gpl" = xno ; then
  AC_MSG_NOTICE([disabled GPL only solvers gpl_solver_methods])
  AC_FOREACH([solver], gpl_solver_methods,[
    if test "x$use_fcs_[]solver[]" = xyes ; then
      use_fcs_[]solver[]="disabled GPL [solver]"
    fi
  ])
fi

#echo "enable_gpl: $enable_gpl"

# add option to disable each solver individually
AC_ARG_ENABLE([fcs-SOLVER],
 [AS_HELP_STRING([--disable-fcs-SOLVER],[disable build of SOLVER (]all_solver_methods[)])])

AC_FOREACH([solver], all_solver_methods,[
  AC_ARG_ENABLE([fcs-]solver,[],,[enable_fcs_]solver[=yes])
  if test "x$enable_fcs_[]solver[]" = xno -a "x$use_fcs_[]solver[]" = xyes ; then
    AC_MSG_NOTICE([[disabling solver method ]]solver)
    use_fcs_[]solver[]="disabled on request"
  fi
])

#AC_FOREACH([solver], all_solver_methods,[
#  echo "[]solver[]: $use_fcs_[]solver[]"
#])

enable_fcs_gridsort=yes # always enabled, because scafacos_test uses gridsort

if test "x$enable_dist" = xyes ; then
  enable_fcs_near=yes
  enable_fcs_resort=yes
fi

# Check prerequisites.

# Common methods and solvers that use these methods.
m4_define([all_common_methods], [fftw3 fftw3_mpi pfft pnfft near gridsort resort])
requires_fftw3="p2nfft p3m"
requires_fftw3_mpi="p2nfft"
requires_pfft="p2nfft"
requires_pnfft="p2nfft"
requires_near="direct ewald p2nfft p3m wolf"
requires_gridsort="ewald memd mmm1d mmm2d p2nfft p3m near"
requires_resort="gridsort"

AC_FOREACH([prereq], all_common_methods,[
  if test "x$enable_fcs_[]prereq[]" = xyes ; then
    use_fcs_[]prereq[]=yes
  else
    use_fcs_[]prereq[]=
    for s in $requires_[]prereq[] ; do
      eval useit=\$use_fcs_$s
      if test "x$useit" = xyes ; then
        use_fcs_[]prereq[]=yes
      fi
    done
  fi
])

if test "x$fcs_float" = "xfloat" ; then
  fft_precision="f"
elif test "x$fcs_float" = "xlong double" ; then
  fft_precision="l"
else
  fft_precision=
fi

# Set up fftw3
use_fcs_fftw3_common=no
use_fcs_fftw3_mpi_common=no
if test "x$use_fcs_fftw3" = xyes -o "x$use_fcs_fftw3_mpi" = xyes ; then
  # Search for installed FFTW3
  if test "x$fcs_float" = "xfloat" ; then
    AX_LIB_FFTW3F_EXT([nobuilt])
  elif test "x$fcs_float" = "xlong double" ; then
    AX_LIB_FFTW3L_EXT([nobuilt])
  else
    AX_LIB_FFTW3_EXT([nobuilt])
  fi

#  echo "use_fcs_fftw3: $use_fcs_fftw3"
#  echo "ax_lib_fftw3: $ax_lib_fftw3"
#  echo "use_fcs_fftw3_mpi: $use_fcs_fftw3_mpi"
#  echo "ax_lib_fftw3_mpi: $ax_lib_fftw3_mpi"

  if test "x$use_fcs_fftw3" = xyes -a "x$ax_lib_fftw3" != xyes || test "x$use_fcs_fftw3_mpi" = xyes -a "x$ax_lib_fftw3_mpi" != xyes || test "x$enable_dist" = xyes ; then

    use_fcs_fftw3_common="$use_fcs_fftw3"
    use_fcs_fftw3_mpi_common="$use_fcs_fftw3_mpi"

    if test "x$use_fcs_fftw3_mpi" = xyes ; then
      AC_MSG_NOTICE([enabling helper method fftw3 with MPI])
      use_fcs_fftw3_mpi="yes (build as helper method)"
    else
      AC_MSG_NOTICE([enabling helper method fftw3])
    fi
    use_fcs_fftw3="yes (build as helper method)"

    # build FFTW without Fortran interface
    case $ac_configure_args in
      *--disable-fortran*) ;;
      *) ac_configure_args="$ac_configure_args --disable-fortran" ;;
    esac

    if test "x$use_fcs_fftw3_mpi_common" = xyes ; then
      # build FFTW with MPI support
      case $ac_configure_args in
        *--enable-mpi*) ;;
        *) ac_configure_args="$ac_configure_args --enable-mpi" ;;
      esac
    fi

    if test "x$fcs_float" = "xfloat" ; then
      # build FFTW for float
      case $ac_configure_args in
        *--enable-float*) ;;
        *) ac_configure_args="$ac_configure_args --enable-float" ;;
      esac
    fi
    if test "x$fcs_float" = "xlong double" ; then
      # build FFTW for long double
      case $ac_configure_args in
        *--enable-long-double*) ;;
        *) ac_configure_args="$ac_configure_args --enable-long-double" ;;
      esac
    fi

    AC_CONFIG_SUBDIRS([lib/common/fftw-3.3])

    AC_DEFINE([FCS_USE_COMMON_FFTW], [1], [Define if solvers should use the FCS internal FFTW3 library.])

    if test -n "$(set `cd "$srcdir" && pwd`; test $# -gt 1 && echo X)" || test -n "$(set `pwd`; test $# -gt 1 && echo X)" ; then
      AC_MSG_ERROR([cannot configure helper method FFTW, because the current directory name and/or the FCS source directory name contains whitespaces!])
    fi

    # Re-check in-tree FFTW3.
    if test "x$fcs_float" = "xfloat" ; then
      AX_LIB_FFTW3F_CHECK([built],[fcs_],[`cd $srcdir && pwd`/lib/common/fftw-3.3/api `cd $srcdir && pwd`/lib/common/fftw-3.3/mpi],[`pwd`/lib/common/fftw-3.3/.libs `pwd`/lib/common/fftw-3.3/mpi/.libs])
    elif test "x$fcs_float" = "xlong double" ; then
      AX_LIB_FFTW3L_CHECK([built],[fcs_],[`cd $srcdir && pwd`/lib/common/fftw-3.3/api `cd $srcdir && pwd`/lib/common/fftw-3.3/mpi],[`pwd`/lib/common/fftw-3.3/.libs `pwd`/lib/common/fftw-3.3/mpi/.libs])
    else
      AX_LIB_FFTW3_CHECK([built],[fcs_],[`cd $srcdir && pwd`/lib/common/fftw-3.3/api `cd $srcdir && pwd`/lib/common/fftw-3.3/mpi],[`pwd`/lib/common/fftw-3.3/.libs `pwd`/lib/common/fftw-3.3/mpi/.libs])
    fi

    # Add in-tree FFTW3 parameters for sub-configures (e.g., pfft, pnfft, p2nfft, ...)
    case $ac_configure_args in
      *--with-fftw3*) ;;
      *)
        ac_configure_args="$ac_configure_args --with-fftw3=built"
        ac_configure_args="$ac_configure_args --with-fftw3-prefix=fcs_"
        ac_configure_args="$ac_configure_args --with-fftw3-includedir=\"`cd $srcdir && pwd`/lib/common/fftw-3.3/api `cd $srcdir && pwd`/lib/common/fftw-3.3/mpi\""
        ac_configure_args="$ac_configure_args --with-fftw3-libdir=\"`pwd`/lib/common/fftw-3.3/.libs `pwd`/lib/common/fftw-3.3/mpi/.libs\""
        ;;
    esac
  fi
fi
AM_CONDITIONAL([ENABLE_COMMON_FFTW],[test "x$use_fcs_fftw3_common" = xyes])

# Set up pfft.
use_fcs_pfft_common=no
if test "x$use_fcs_pfft" = xyes ; then
  # Search for installed PFFT.
  if test "x$fcs_float" = "xfloat" ; then
    AX_LIB_PFFTF_EXT([nobuilt])
  elif test "x$fcs_float" = "xlong double" ; then
    AX_LIB_PFFTL_EXT([nobuilt])
  else
    AX_LIB_PFFT_EXT([nobuilt])
  fi

#  echo "use_fcs_pfft: $use_fcs_pfft"
#  echo "ax_lib_pfft: $ax_lib_pfft"

  if test "x$ax_lib_pfft" != xyes || test "x$enable_dist" = xyes ; then
    use_fcs_pfft_common=yes

    AC_MSG_NOTICE([enabling helper method pfft])

    use_fcs_pfft="yes (build as helper method)"

    # Do not build shared libraries, ScaFaCoS is not ready for it.
    case $ac_configure_args in
      *--disable-shared*) ;;
      *) ac_configure_args="$ac_configure_args --disable-shared" ;;
    esac
    # Disable install of headers.
    case $ac_configure_args in
      *--disable-header-install*) ;;
      *) ac_configure_args="$ac_configure_args --disable-header-install" ;;
    esac
    # Disable compilation of documentation.
    case $ac_configure_args in
      *--disable-doc*) ;;
      *) ac_configure_args="$ac_configure_args --disable-doc" ;;
    esac
    # Disable tests if extended tests are not requested.
    if test "x$enable_fcs_extended_tests" != xyes ; then
      case $ac_configure_args in
        *--disable-tests*) ;;
        *) ac_configure_args="$ac_configure_args --disable-tests" ;;
      esac
    fi

    if test "x$fcs_float" = "xfloat" ; then
      # build PFFT for float
      case $ac_configure_args in
        *--enable-float*) ;;
        *) ac_configure_args="$ac_configure_args --enable-float" ;;
      esac
    fi
    if test "x$fcs_float" = "xlong double" ; then
      # build PFFT for long double
      case $ac_configure_args in
        *--enable-long-double*) ;;
        *) ac_configure_args="$ac_configure_args --enable-long-double" ;;
      esac
    fi

    # Enable fcs_ prefix.
    case $ac_configure_args in
      *--enable-lib-prefix*) ;;
      *) ac_configure_args="$ac_configure_args --enable-lib-prefix=fcs_" ;;
    esac

    AC_CONFIG_SUBDIRS([lib/common/pfft])

    AC_DEFINE([FCS_USE_COMMON_PFFT], [1], [Define if solvers should use the FCS internal PFFT library.])

    if test -n "$(set `cd "$srcdir" && pwd`; test $# -gt 1 && echo X)" || test -n "$(set `pwd`; test $# -gt 1 && echo X)" ; then
      AC_MSG_ERROR([cannot configure helper method PFFT, because the current directory name and/or the FCS source directory name contains whitespaces!])
    fi

    # Re-check in-tree PFFT.
    if test "x$fcs_float" = "xfloat" ; then
      AX_LIB_PFFTF_CHECK([built],[],[`cd $srcdir && pwd`/lib/common/pfft/api],[`pwd`/lib/common/pfft/.libs])
    elif test "x$fcs_float" = "xlong double" ; then
      AX_LIB_PFFTL_CHECK([built],[],[`cd $srcdir && pwd`/lib/common/pfft/api],[`pwd`/lib/common/pfft/.libs])
    else
      AX_LIB_PFFT_CHECK([built],[],[`cd $srcdir && pwd`/lib/common/pfft/api],[`pwd`/lib/common/pfft/.libs])
    fi

    # Add in-tree PFFT parameters for sub-configures (e.g., pnfft, p2nfft, ...)
    case $ac_configure_args in
      *--with-pfft*) ;;
      *)
        ac_configure_args="$ac_configure_args --with-pfft=built"
        ac_configure_args="$ac_configure_args --with-pfft-prefix=fcs_"
        ac_configure_args="$ac_configure_args --with-pfft-includedir=\"`cd $srcdir && pwd`/lib/common/pfft/api\""
        ac_configure_args="$ac_configure_args --with-pfft-libdir=\"`pwd`/lib/common/pfft/.libs\""
        ;;
    esac
  fi
fi
AM_CONDITIONAL([ENABLE_COMMON_PFFT],[test "x$use_fcs_pfft_common" = xyes])

# Set up pnfft.
use_fcs_pnfft_common=no
if test "x$use_fcs_pnfft" = xyes ; then
  # Search for installed PNFFT.
  if test "x$fcs_float" = "xfloat" ; then
    AX_LIB_PNFFTF_EXT([nobuilt])
  elif test "x$fcs_float" = "xlong double" ; then
    AX_LIB_PNFFTL_EXT([nobuilt])
  else
    AX_LIB_PNFFT_EXT([nobuilt])
  fi

#  echo "use_fcs_pnfft: $use_fcs_pnfft"
#  echo "ax_lib_pnfft: $ax_lib_pnfft"

  if test "x$ax_lib_pnfft" != xyes || test "x$enable_dist" = xyes ; then
    use_fcs_pnfft_common=yes

    AC_MSG_NOTICE([enabling helper method pnfft])

    use_fcs_pnfft="yes (build as helper method)"

    # Do not build shared libraries, ScaFaCoS is not ready for it.
    case $ac_configure_args in
      *--disable-shared*) ;;
      *) ac_configure_args="$ac_configure_args --disable-shared" ;;
    esac
    # Disable install of headers.
    case $ac_configure_args in
      *--disable-header-install*) ;;
      *) ac_configure_args="$ac_configure_args --disable-header-install" ;;
    esac
    # Disable compilation of documentation.
    case $ac_configure_args in
      *--disable-doc*) ;;
      *) ac_configure_args="$ac_configure_args --disable-doc" ;;
    esac
    # Disable tests if extended tests are not requested.
    if test "x$enable_fcs_extended_tests" != xyes ; then
      case $ac_configure_args in
        *--disable-tests*) ;;
        *) ac_configure_args="$ac_configure_args --disable-tests" ;;
      esac
    fi

    if test "x$fcs_float" = "xfloat" ; then
      # build PNFFT for float
      case $ac_configure_args in
        *--enable-float*) ;;
        *) ac_configure_args="$ac_configure_args --enable-float" ;;
      esac
    fi
    if test "x$fcs_float" = "xlong double" ; then
      # build PNFFT for long double
      case $ac_configure_args in
        *--enable-long-double*) ;;
        *) ac_configure_args="$ac_configure_args --enable-long-double" ;;
      esac
    fi

    # Enable fcs_ prefix.
    case $ac_configure_args in
      *--enable-lib-prefix*) ;;
      *) ac_configure_args="$ac_configure_args --enable-lib-prefix=fcs_" ;;
    esac

    AC_CONFIG_SUBDIRS([lib/common/pnfft])

    AC_DEFINE([FCS_USE_COMMON_PNFFT], [1], [Define if solvers should use the FCS internal PNFFT library.])

    if test -n "$(set `cd "$srcdir" && pwd`; test $# -gt 1 && echo X)" || test -n "$(set `pwd`; test $# -gt 1 && echo X)" ; then
      AC_MSG_ERROR([cannot configure helper method PNFFT, because the current directory name and/or the FCS source directory name contains whitespaces!])
    fi

    # Re-check in-tree PNFFT.
    if test "x$fcs_float" = "xfloat" ; then
      AX_LIB_PNFFTF_CHECK([built],[],[`cd $srcdir && pwd`/lib/common/pnfft/api],[`pwd`/lib/common/pnfft/.libs])
    elif test "x$fcs_float" = "xlong double" ; then
      AX_LIB_PNFFTL_CHECK([built],[],[`cd $srcdir && pwd`/lib/common/pnfft/api],[`pwd`/lib/common/pnfft/.libs])
    else
      AX_LIB_PNFFT_CHECK([built],[],[`cd $srcdir && pwd`/lib/common/pnfft/api],[`pwd`/lib/common/pnfft/.libs])
    fi

    # Add in-tree PNFFT parameters for sub-configures (e.g., p2nfft, ...)
    case $ac_configure_args in
      *--with-pnfft*) ;;
      *)
        ac_configure_args="$ac_configure_args --with-pnfft=built"
        ac_configure_args="$ac_configure_args --with-pnfft-prefix=fcs_"
        ac_configure_args="$ac_configure_args --with-pnfft-includedir=\"`cd $srcdir && pwd`/lib/common/pnfft/api\""
        ac_configure_args="$ac_configure_args --with-pnfft-libdir=\"`pwd`/lib/common/pnfft/.libs\""
        ;;
    esac
  fi
fi
AM_CONDITIONAL([ENABLE_COMMON_PNFFT],[test "x$use_fcs_pnfft_common" = xyes])

# Set up near field module.
if test "x$use_fcs_near" = xyes ; then
  AC_MSG_NOTICE([enabling helper method 'near'])
fi
AM_CONDITIONAL([ENABLE_COMMON_NEAR],[test "x$use_fcs_near" = xyes])

# Set up gridsort module.
if test "x$use_fcs_gridsort" = xyes ; then
  AC_MSG_NOTICE([enabling helper method 'gridsort'])
fi
AM_CONDITIONAL([ENABLE_COMMON_GRIDSORT],[test "x$use_fcs_gridsort" = xyes])

# Set up resort module.
if test "x$use_fcs_resort" = xyes ; then
  AC_MSG_NOTICE([enabling helper method 'resort'])
fi
AM_CONDITIONAL([ENABLE_COMMON_RESORT],[test "x$use_fcs_resort" = xyes])

# Solver prerequisites.

# direct prerequisites
AX_FCS_DIRECT_TOP([test "x$use_fcs_direct" = xyes])

# ewald prerequisites

# fmm prerequisites
AX_FCS_FMM_TOP([test "x$use_fcs_fmm" = xyes])

# p2nfft prerequisites

# p3m prerequisites

# pepc prerequisites

# pp3mg prerequisites

# vmg prerequisites

# wolf prerequisites

enabled_solvers=

# add Automake conditional ENABLE_<SOLVER> for each SOLVER
AC_FOREACH([solver], all_solver_methods,[
  AM_CONDITIONAL([ENABLE_]m4_toupper(solver), [test "x$use_fcs_]solver[" = xyes])
  if test "x$use_fcs_[]solver" = xyes ; then
    AC_DEFINE_UNQUOTED([FCS_ENABLE_]m4_toupper(solver), [1],
      [Whether solver method ]solver[ is enabled.])
    enabled_solvers="${enabled_solvers},solver"
  fi
])

enabled_solvers="${enabled_solvers#,}"

##################################################
# SUPPORT META-PACKAGE WITH REQUIRED INFORMATION.
# common methods
AX_FCS_PACKAGE_RESET
AX_FCS_PACKAGE_ADD([fcs_USE],[yes])
if test "x$use_fcs_fortran" = xyes ; then
  AX_FCS_PACKAGE_ADD([fcs_LIBS],[-lfcs4fortran])
  AX_FCS_PACKAGE_ADD([fcs_LIBS_A],[src/fortran/libfcs4fortran.a])
fi
AX_FCS_PACKAGE_ADD([fcs_LIBS],[-lfcs])
AX_FCS_PACKAGE_ADD([fcs_LIBS_A],[src/libfcs.a])
if test "x$use_fcs_near" = xyes ; then
  AX_FCS_PACKAGE_ADD([near_USE],[yes])
  AX_FCS_PACKAGE_ADD([near_LIBS],[-lfcs_near])
  AX_FCS_PACKAGE_ADD([near_LIBS_A],[lib/common/near/libfcs_near.a])
fi
if test "x$use_fcs_gridsort" = xyes ; then
  AX_FCS_PACKAGE_ADD([gridsort_USE],[yes])
  AX_FCS_PACKAGE_ADD([gridsort_LIBS],[-lfcs_gridsort])
  AX_FCS_PACKAGE_ADD([gridsort_LIBS_A],[lib/common/gridsort/libfcs_gridsort.a])
fi
if test "x$use_fcs_resort" = xyes ; then
  AX_FCS_PACKAGE_ADD([resort_USE],[yes])
  AX_FCS_PACKAGE_ADD([resort_LIBS],[-lfcs_resort])
  AX_FCS_PACKAGE_ADD([resort_LIBS_A],[lib/common/resort/libfcs_resort.a])
fi
AX_FCS_PACKAGE_ADD([fcs_common_USE],[yes])
AX_FCS_PACKAGE_ADD([fcs_common_LIBS],[-lfcs_common])
AX_FCS_PACKAGE_ADD([fcs_common_LIBS_A],[lib/common/fcs-common/libfcs_common.a])
if test "x$use_fcs_pfft_common" = xyes ; then
  AX_FCS_PACKAGE_ADD([pfft_USE],[yes])
  AX_FCS_PACKAGE_ADD([pfft_LIBS],[-lfcs_pfft${fft_precision}])
  AX_FCS_PACKAGE_ADD([pfft_LIBS_A],[lib/common/pfft/.libs/libfcs_pfft${fft_precision}.a])
fi
if test "x$use_fcs_pnfft_common" = xyes ; then
  AX_FCS_PACKAGE_ADD([pnfft_USE],[yes])
  AX_FCS_PACKAGE_ADD([pnfft_LIBS],[-lfcs_pnfft${fft_precision}])
  AX_FCS_PACKAGE_ADD([pnfft_LIBS_A],[lib/common/pnfft/.libs/libfcs_pnfft${fft_precision}.a])
fi
if test "x$use_fcs_fftw3_common" = xyes ; then
  AX_FCS_PACKAGE_ADD([fftw3_common_USE],[yes])
  AX_FCS_PACKAGE_ADD([fftw3_common_LIBS],[-lfcs_fftw3${fft_precision}])
  AX_FCS_PACKAGE_ADD([fftw3_common_LIBS_A],[lib/common/fftw-3.3/.libs/libfcs_fftw3${fft_precision}.a])
fi
if test "x$use_fcs_fftw3_mpi_common" = xyes ; then
  AX_FCS_PACKAGE_ADD([fftw3_mpi_common_USE],[yes])
  AX_FCS_PACKAGE_ADD([fftw3_mpi_common_LIBS],[-lfcs_fftw3${fft_precision}_mpi])
  AX_FCS_PACKAGE_ADD([fftw3_mpi_common_LIBS_A],[lib/common/fftw-3.3/mpi/.libs/libfcs_fftw3${fft_precision}_mpi.a])
fi
if test "x$use_fcs_fftw3" = xyes ; then
  AX_FCS_PACKAGE_ADD([fftw3_USE],[yes])
  AX_FCS_PACKAGE_ADD([fftw3_LDADD],[${fftw3_LDFLAGS} ${fftw3_LIBS}])
fi
if test "x$use_fcs_fftw3_mpi" = xyes ; then
  AX_FCS_PACKAGE_ADD([fftw3_mpi_USE],[yes])
  AX_FCS_PACKAGE_ADD([fftw3_mpi_LDADD],[${fftw3_LDFLAGS} ${fftw3_mpi_LIBS}])
fi


# specific solvers
if test "x$use_fcs_direct" = xyes ; then
  AC_CONFIG_FILES([lib/direct/Makefile])
  AX_FCS_PACKAGE_ADD([direct_LIBS],[-lfcs_direct])
  AX_FCS_PACKAGE_ADD([direct_LIBS_A],[lib/direct/libfcs_direct.a])
fi
if test "x$use_fcs_memd" = xyes ; then
  AC_CONFIG_FILES([lib/memd/Makefile])
  AX_FCS_PACKAGE_ADD([memd_LIBS],[-lfcs_memd])
  AX_FCS_PACKAGE_ADD([memd_LIBS_A],[lib/memd/libfcs_memd.a])
fi
if test "x$use_fcs_mmm1d" = xyes ; then
  AC_CONFIG_FILES([lib/mmm1d/Makefile])
  AX_FCS_PACKAGE_ADD([mmm1d_LIBS],[-lfcs_mmm1d])
  AX_FCS_PACKAGE_ADD([mmm1d_LIBS_A],[lib/mmm1d/libfcs_mmm1d.a])
fi
if test "x$use_fcs_mmm2d" = xyes ; then
  AC_CONFIG_FILES([lib/mmm2d/Makefile])
  AX_FCS_PACKAGE_ADD([mmm2d_LIBS],[-lfcs_mmm2d])
  AX_FCS_PACKAGE_ADD([mmm2d_LIBS_A],[lib/mmm2d/libfcs_mmm2d.a])
fi
if test "x$use_fcs_p3m" = xyes ; then
  AC_CONFIG_FILES([lib/p3m/Makefile])
  AX_FCS_PACKAGE_ADD([p3m_LIBS],[-lfcs_p3m])
  AX_FCS_PACKAGE_ADD([p3m_LIBS_A],[lib/p3m/libfcs_p3m.a])
  AX_FCS_PACKAGE_ADD([fftw3_USE],[yes])
fi
if test "x$use_fcs_ewald" = xyes ; then
  AC_CONFIG_FILES([lib/ewald/Makefile])
  AX_FCS_PACKAGE_ADD([ewald_LIBS],[-lfcs_ewald])
  AX_FCS_PACKAGE_ADD([ewald_LIBS_A],[lib/ewald/libfcs_ewald.a])
fi
if test "x$use_fcs_wolf" = xyes ; then
  AC_CONFIG_FILES([lib/wolf/Makefile])
  AX_FCS_PACKAGE_ADD([wolf_LIBS],[-lfcs_wolf])
  AX_FCS_PACKAGE_ADD([wolf_LIBS_A],[lib/wolf/libfcs_wolf.a])
fi

# Create FCS library files.
AC_CONFIG_FILES([Makefile
  lib/Makefile
  lib/common/fcs-common/Makefile
  src/Makefile
  src/fortran/Makefile])

# Create near-field module files.
AC_CONFIG_FILES([lib/common/near/Makefile])

# Create gridsort module files.
AC_CONFIG_FILES([lib/common/gridsort/Makefile])

# Create resort module files.
AC_CONFIG_FILES([lib/common/resort/Makefile])

# Create test files.
AC_CONFIG_FILES([
  test/Makefile 
  test/generic/Makefile 
  test/c/Makefile 
  test/fortran/Makefile
])

AC_CONFIG_LINKS([test/generic/scafacos_test_all.sh:test/generic/scafacos_test_all.sh])
AC_CONFIG_LINKS([test/generic/compute_nonperiodic_references.sh:test/generic/compute_nonperiodic_references.sh])
AC_CONFIG_LINKS([test/generic/opt_alpha_cloud_wall_2dp.py:test/generic/opt_alpha_cloud_wall_2dp.py])
AC_CONFIG_LINKS([test/generic/opt_alpha_cloud_wall_3dp.py:test/generic/opt_alpha_cloud_wall_3dp.py])

# Create helper script for portable MPI job startup.
AC_CONFIG_FILES([test/defs])

##################################################
# DOCUMENTATION
#
AC_ARG_ENABLE([doc],
 [AS_HELP_STRING([--disable-doc],[disable build of documents (PDF manual etc.)])],,
 [enable_doc=yes])

if test "x$enable_doc" != xno -a ! -d "$srcdir/doc"; then
  AC_MSG_WARN([disabling doc as it is not present in source tree])
  enable_doc=no
fi

if test "x$enable_doc" != xno ; then
  AC_PATH_PROG(PDFLATEX,pdflatex,no)
  AC_PATH_PROG(MAKEINDEX,makeindex,no)
  AC_PATH_PROG(BIBTEX,bibtex,no)
  AC_CHECK_PROGS([DOXYGEN], [doxygen])
    if test -z "$DOXYGEN";
       then AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
    fi
else
  AC_MSG_NOTICE([disable build of documents])
fi

AM_CONDITIONAL(ENABLE_DOC,[test "x$enable_doc" != xno])
AM_CONDITIONAL(HAVE_LATEX,
	[test "x$PDFLATEX" != xno && 
	test "x$MAKEINDEX" != xno && 
	test "x$BIBTEX" != xno])

if test "x$enable_doc" != xno ; then
  AC_CONFIG_FILES([doc/Makefile doc/latexinc.sh])
  AC_CONFIG_LINKS([doc/latex-manual.sh:doc/latex-manual.sh])
fi

AM_CONDITIONAL([HAVE_DOXYGEN], [test -n "$DOXYGEN"])
AM_COND_IF([HAVE_DOXYGEN], [AC_CONFIG_FILES([doc/Doxyfile])])

# Create config headers: config.h comes from autoheader-generated config.h.in,
# fcs_config.h.in is written manually and contains the public parts of config.h,
# fconfig.h.in is written manually and only contains stuff needed for Fortran.
AC_CONFIG_HEADERS([config.h fcs_config.h])
# Convert fconfig.h.in into an fconfig.h suitable for inclusion into Fortran.
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
# Convert fcs_fconfig.h.in into fcs_fconfig.h suitable for inclusion into Fortran.
# The sed script transforms C style /* comments */ into Fortran style.
AC_CONFIG_HEADERS([fcs_fconfig.h],
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
    }' fcs_fconfig.h > fcs_fconfig.ht
   mv -f fcs_fconfig.ht fcs_fconfig.h
]])

# Descend into sub packages.
AC_FOREACH([solver], subconfigure_solver_methods,
  [if test "x$use_fcs_[]solver[]" = xyes ; then
    AC_CONFIG_SUBDIRS([lib/]solver)
   fi
])

# Decend into meta package that gathers informations from all other packages.
ac_configure_args="$ac_configure_args --enable-fcs-package-solvers=${enabled_solvers}"
AC_CONFIG_SUBDIRS([package])

AC_OUTPUT

AC_MSG_NOTICE([])
AC_MSG_NOTICE([****************************************************************])
AC_MSG_NOTICE([*                   Configuration summary                      *])
AC_MSG_NOTICE([****************************************************************])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([Compiler settings:])
AC_MSG_NOTICE([------------------])
AC_MSG_NOTICE([ C compiler                 CC        = $CC])
AC_MSG_NOTICE([ C compiler flags           CFLAGS    = $CFLAGS])
AC_MSG_NOTICE([ C preprocessor flags       CPPFLAGS  = $CPPFLAGS])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([ C++ compiler               CXX       = $CXX])
AC_MSG_NOTICE([ C++ compiler flags         CXXFLAGS  = $CXXFLAGS])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([ Fortran compiler           FC        = $FC])
AC_MSG_NOTICE([ Fortran compiler flags     FCFLAGS   = $FCFLAGS])
AC_MSG_NOTICE([ Fortran preprocessor flags FCPPFLAGS = $FCPPFLAGS])
AC_MSG_NOTICE([ Fortran definitions        FCDEFS    = $FCDEFS])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([ Libraries                  LIBS      = $LIBS])
AC_MSG_NOTICE([])
#AC_MSG_NOTICE([Requirements:])
#AC_MSG_NOTICE([-------------])
#AC_MSG_NOTICE([])
AC_MSG_NOTICE([Package settings:])
AC_MSG_NOTICE([-----------------])
AC_MSG_NOTICE([ FCS C types                = $fcs_int / $fcs_float])
AC_MSG_NOTICE([ FCS Fortan types           = $fcs_integer / $fcs_real])
AC_MSG_NOTICE([ Fortan interface           = $use_fcs_fortran])
AC_MSG_NOTICE([])
AC_FOREACH([method], all_common_methods,
[
m=`printf %26s method`
m="method${m%method}"
if test "x${use_fcs_[]method[]%% *}" = xyes ; then
  AC_MSG_NOTICE([ ${m} = $use_fcs_[]method[]])
elif test -n "$use_fcs_[]method[]" ; then
  AC_MSG_NOTICE([ ${m} = no ($use_fcs_[]method)])
else
  AC_MSG_NOTICE([ ${m} = no])
fi
])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([ Use GPL solvers            = $enable_gpl])
AC_FOREACH([method], all_solver_methods,
[
m=`printf %26s method`
m="method${m%method}"
if test "x$use_fcs_[]method[]" = xyes ; then
  AC_MSG_NOTICE([ ${m} = yes])
elif test -n "$use_fcs_[]method[]" ; then
  AC_MSG_NOTICE([ ${m} = no ($use_fcs_[]method)])
else
  AC_MSG_NOTICE([ ${m} = no])
fi
])

# Check whether there is at least one solver available!
if AC_FOREACH([solver], all_solver_methods, [test "x$use_fcs_[]solver[]" != xyes && ]) true ; then
  AC_MSG_WARN([no valid solver specified, or missing prerequisites])
  AC_MSG_ERROR([need to specify at least one valid solver.
		Available solvers are: ]all_solver_methods)
fi
