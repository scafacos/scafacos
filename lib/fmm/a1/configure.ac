dnl
dnl $Id: configure.in 16 2008-02-05 02:18:30Z balaji $
dnl
dnl Copyright (C) 2010 by Argonne National Laboratory.
dnl     See COPYRIGHT in top-level directory.
dnl

AC_PREREQ(2.62)

AC_INIT([a1],[0])

AC_MSG_NOTICE([****************************************************************])
AC_MSG_NOTICE([*      Configuring A1 (part of solver FMM)                     *])
AC_MSG_NOTICE([****************************************************************])

AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR(m4)
AM_INIT_AUTOMAKE([-Wall -Werror foreign 1.11 silent-rules subdir-objects])

LT_PREREQ([2.2.6])

AX_PROG_CC_MPI(,,AC_MSG_FAILURE([The A1 library requires an MPI C compiler.]))

# Bug in libtool adds -O2 and -g by default
# uncomment for non-POSIX systems
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
PAC_PUSH_FLAG(CFLAGS)
LT_INIT(disable-shared)
PAC_POP_FLAG(CFLAGS)

if test ! -s "$srcdir/maint/version" ; then
   AC_ERROR([Version information not found])
fi

# ABI version
liba1_abi_version=`cat $srcdir/maint/version | grep ^liba1_abi_version | cut -f3 -d' '`
AC_SUBST(liba1_abi_version)

# Release version
# Produce a numeric version assuming the following format:
# Version: [MAJ].[MIN].[REV][EXT][EXT_NUMBER]
# Example: 1.0.7rc1 has
#          MAJ = 1
#          MIN = 0
#          REV = 7
#          EXT = rc
#          EXT_NUMBER = 1
#
# Converting to numeric version will convert EXT to a format number:
#          ALPHA (a) = 0
#          BETA (b)  = 1
#          RC (rc)   = 2
#          PATCH (p) = 3
# Regular releases are treated as patch 0
#
# Numeric version will have 1 digit for MAJ, 2 digits for MIN,
# 2 digits for REV, 1 digit for EXT and 2 digits for EXT_NUMBER.
A1_VERSION=`cat $srcdir/maint/version | grep ^A1_VERSION | cut -f3 -d' '`
changequote(<<,>>)
V1=`expr $A1_VERSION : '\([0-9]*\)\.[0-9]*\.*[0-9]*[a-zA-Z]*[0-9]*'`
V2=`expr $A1_VERSION : '[0-9]*\.\([0-9]*\)\.*[0-9]*[a-zA-Z]*[0-9]*'`
V3=`expr $A1_VERSION : '[0-9]*\.[0-9]*\.*\([0-9]*\)[a-zA-Z]*[0-9]*'`
V4=`expr $A1_VERSION : '[0-9]*\.[0-9]*\.*[0-9]*\([a-zA-Z]*\)[0-9]*'`
V5=`expr $A1_VERSION : '[0-9]*\.[0-9]*\.*[0-9]*[a-zA-Z]*\([0-9]*\)'`
changequote([,])

if test "$V2" -le 9 ; then V2=0$V2 ; fi
if test "$V3" = "" ; then V3=0; fi
if test "$V3" -le 9 ; then V3=0$V3 ; fi
if test "$V4" = "a" ; then
    V4=0
elif test "$V4" = "b" ; then
    V4=1
elif test "$V4" = "rc" ; then
    V4=2
elif test "$V4" = "" ; then
    V4=3
    V5=0
elif test "$V4" = "p" ; then
    V4=3
fi
if test "$V5" -le 9 ; then V5=0$V5 ; fi

A1_NUMVERSION=`expr $V1$V2$V3$V4$V5 + 0`
AC_SUBST(A1_NUMVERSION)
AC_SUBST(A1_VERSION)

AC_CONFIG_HEADER(src/include/a1conf.h)
AH_TOP([/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  (C) 2010 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */
#ifndef A1CONF_H_INCLUDED
#define A1CONF_H_INCLUDED
])
AH_BOTTOM([#endif])

## Some of the autoconf macros automatically edit the CFLAGS and
## friends. We save and reset them to avoid this behavior. We want to
## control what flags get set.
#PAC_PUSH_ALL_FLAGS
#AC_PROG_CC
#PAC_POP_ALL_FLAGS

PAC_ARG_STRICT
PAC_CC_FUNCTION_NAME_SYMBOL

## Check if the necessary headers are available
AC_CHECK_HEADERS(stdio.h stdlib.h string.h strings.h unistd.h stdarg.h sys/types.h \
		sys/stat.h time.h errno.h)

## Check what device we want to use
a1_device_dcmfd=false
a1_device_pthreadd=false
a1_device_mpid=false
a1_device_geminid=false
a1_device_lapid=false
a1_device_ibd=false
#AC_ARG_WITH(device, AC_HELP_STRING([--with-device],[A1 device to use]),
#	a1_device=$withval,a1_device=dcmfd)
a1_device=dcmfd
got_device=no
for dev in `echo $a1_device | sed -e 's/:/ /g' -e 's/,/ /g'` ; do
    case "$dev" in
    	 dcmfd)
		PAC_SET_HEADER_LIB_PATH(dcmf)
                PAC_SET_HEADER_LIB_PATH(spi)
		PAC_CHECK_HEADER_LIB(spi,spi/DMA_Addressing.h,SPI.cna,DMA_AddressingGetHandle,
				got_spi=yes,got_spi=no)
		got_dcmf=no
		if test "$got_spi" = "yes" ; then
			PAC_CHECK_HEADER_LIB(dcmf,dcmf.h,dcmf-fast.cnk,DCMF_Messager_initialize,
				got_dcmf=yes,got_dcmf=no)
		fi
		if test "$got_dcmf" = "yes" ; then
		   got_device=yes
                   a1_device_dcmfd=true
		fi
                got_device=yes
                a1_device_dcmfd=true
		;;
#         pthreadd)
#        PAC_SET_HEADER_LIB_PATH(pthread)
#        got_pthread=no
#        PAC_CHECK_HEADER_LIB(pthread,pthread.h,pthread,pthread_create,
#                got_pthread=yes,got_pthread=no)
#        if test "$got_pthread" = "yes" ; then
#           got_device=yes
#                   a1_device_pthreadd=true
#        fi
#                got_device=yes
#                a1_device_pthreadd=true
#        ;;
#    	 mpid)
#		PAC_SET_HEADER_LIB_PATH(mpi)
#		got_mpi=no
#		PAC_CHECK_HEADER_LIB(mpich,mpi.h,mpich,MPI_Init_thread,
#			got_mpi=yes,got_mpi=no)
#		if test "$got_mpi" = "yes" ; then
#		   got_device=yes
#                   a1_device_mpid=true
#		fi
#                got_device=yes
#                a1_device_mpid=true
#		;;
	 *)
		;;
    esac
done
AM_CONDITIONAL([a1_device_dcmfd],[$a1_device_dcmfd])
if test "$got_device" = "no" ; then
   AC_ERROR([no valid device found])
fi

## Error checking functionality
AC_ARG_ENABLE(error-checking,
	AC_HELP_STRING([--enable-error-checking],[Enable error checking functionality]),
	enable_error_checking=$enableval,
	enable_error_checking=yes)
if test "$enable_error_checking" = "yes" ; then
   AC_DEFINE(ERROR_CHECKING,1,[Define if error checking is enabled])
fi

## Check if __VA_ARGS__ is defined by the compiler
PAC_C_MACRO_VA_ARGS

## Strict compilation
PAC_ARG_STRICT

## Debugging support
AC_ARG_ENABLE(g, AC_HELP_STRING([--enable-g],[Enable Debugging]),
                 [ debug=$enableval ],
                 [ debug=no ])
AC_MSG_CHECKING(debugging support)
AC_MSG_RESULT($debug)
if test "$debug" = "yes"; then
   CFLAGS="$CFLAGS -g"
   AC_DEFINE(ENABLE_DEBUG,1,[Define if debugging support is enabled])
fi

## Disable alignment checking
AC_ARG_ENABLE(alignment-checking,
	AC_HELP_STRING([--enable-alignment-checking],[Enable alignment checking]),
                 [ align_check=$enableval ],
                 [ align_check=yes ])
AC_MSG_CHECKING(alignment checking)
AC_MSG_RESULT($align_check)
if test "$align_check" = "yes"; then
   AC_DEFINE(ENABLE_ALIGNMENT_CHECKING,1,[Define if alignment checking is enabled])
fi

## Documentation
AC_PATH_PROG([DOXYGEN],[doxygen],,$PATH)
AC_SUBST(DOXYGEN)

#AC_SUBST(top_srcdir)

m4_ifdef([AX_FCS_PACKAGE_RESET],[
AX_FCS_PACKAGE_RESET
AX_FCS_PACKAGE_ADD([fmm_LDADD],[$LDFLAGS $LIBS])
])

## Final output
AC_OUTPUT(Makefile
	src/include/a1.h
	a1-doxygen.cfg
	)
