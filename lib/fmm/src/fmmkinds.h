#ifdef HAVE_FCONFIG_H
#include <fconfig.h>
#endif

! define common FMM kind types
#ifdef fcs_real_kind
#define FMM_REAL fcs_real_kind
#define FMM_REAL_EXTENDED fcs_real_kind
#else
#define FMM_REAL 8
#define FMM_REAL_EXTENDED 8
#endif
#define FMM_REAL_ITOR 8

!FMM needs 8-byte integer types internally
#if (fcs_integer_kind > 8)
#define FMM_INTEGER fcs_integer_kind
#define FMM_LOGICAL fcs_integer_kind
#else
#define FMM_INTEGER 8
#define FMM_LOGICAL 8
#endif

#ifdef fcs_integer_kind_isoc
#define FMM_INT_ISOC_KIND fcs_integer_kind_isoc
#else
#define FMM_INT_ISOC_KIND c_int
#endif

#ifdef fcs_real_kind_isoc
#define FMM_REAL_ISOC_KIND fcs_real_kind_isoc
#else
#define FMM_REAL_ISOC_KIND c_double
#endif

! generate correct C types from Fortran types for iso binding
#if FMM_INTEGER == 8
#define FMM_C_INTEGER c_long_long
#elif FMM_INTEGER == 4
#define FMM_C_INTEGER c_int
#else
#error "fmm integer type not supported."
#endif

#if FMM_REAL == 16
#define FMM_C_REAL c_long_double
#elif FMM_REAL == 8
#define FMM_C_REAL c_double
#elif FMM_REAL == 4
#define FMM_C_REAL c_float
#else
#error "fmm real type not supported."
#endif

! define C ptr size for Fortran
!#define FMM_POINTER SIZEOF_INT_P
!TODO 
#define FMM_POINTER 4

!allow bit manipulation of floats
#if (FMM_REAL == 4)
#define FMM_XYZ_TO_INTEGER  4
#elif (FMM_REAL == 8)
#define FMM_XYZ_TO_INTEGER  8
#else
#define FMM_XYZ_TO_INTEGER -1
#endif

! define debug allocation kind type
#define FMM_TESTALLOC_INTEGER 8

! additional statements for the parallel code
!
! maximum number of processes = 2**(8*FMM_MP_INTEGER_PROCESSES-1)
! ======================================================================
! | FMM_MP_INTEGER_PROCESSES |       maximum number of processes       |
! ======================================================================
! |             1            |                   128                   |
! |             2            |                  32768                  |
! |             4            |                2147483648               |
! |             8            |           9223372036854775808           |
! |            16            | 170141183460469231731687303715884105728 |
! ======================================================================

#ifdef fcs_real_kind
#define FMM_MP_REAL_MAX fcs_real_kind
#else
#define FMM_MP_REAL_MAX 8
#endif
#define FMM_MP_INTEGER_PROCESSES 8 

#ifdef FMM_COMPRESSION
#ifdef FMM_PARALLEL
#ifdef FMM_LOADSORT
#error "compression & loadsort not supported."
#endif
#endif
#endif

