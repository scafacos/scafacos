! get the defines computed by configure
#include <fconfig.h>

! limit number of multipoles
!set by Scafacos configure 

! enable parallel FMM
#define FMM_PARALLEL

! print output to console
#undef FMM_INFO

! enable debug mode
#undef FMM_DEBUG

! enable time measurement, walltime OR cputime
#define FMM_CPUTIME
#undef FMM_WALLTIME

! print additional FMM information
#define FMM_STATISTICS

! enable improved charge-neutral error estimation
#define FMM_DAMPING
! compute ++ energy far field in dipole approximation
#define FMM_DAMPING_PP_DIPOLEMOMENTS

! diable potential vector
#undef FMM_NOPOT

! disable gradient vector
#define FMM_NOGRAD

! correct center of mass movement
#define FMM_CORRECTION_OF_FORCES

! Restore input coordinates after FMM run
#define FMM_RESTORE_COORDINATES

! allow allocation of additional tree management buffer
#define FMM_IBOXSCR

! allow additional FMM-defined memory for libsl
#define FMM_SORTMEMORY

! enable alternative sequential sorting algorithm
#undef FMM_SORTHD
#undef FMM_SORTHDM

! enable tree-like summation of multipole moments (esp. for single precision)
#undef FMM_MULTIPOLEMOMENTS

! enable fortran allocation test routines
! for local allocation
#undef FMM_ALLOCATION
! for global allocation
#undef FMM_ALLOCATIONALL

! enable aligned memory allocation
#define FMM_ALLOCALIGNED
#define FMM_ALIGNMENT 16

! disable F2003 function pointer features
!#undef FMM_NOFUNCTIONPOINTER ! this is now set by configure (M. Hofmann)

! enable c bindings interoperability feature
#define FMM_ISO_C_BINDING

!enable constants from C source files for reduced compile time
#define FMM_C_CONSTANTS
#undef FMM_C_CONSTANTS_TODO

! enable memory compression scheme (esp. for single precision, low precision)
#undef FMM_TREETOGRAD
#undef FMM_EXTREMETREETOGRAD
#undef FMM_COMPRESSION
#undef FMM_EXTREMECOMPRESSION
#undef FMM_EXTREMEEXTREMECOMPRESSION
#undef FMM_SIGNEXPONENT

! enable parallel load balancing (nyi!)
#undef FMM_LOADSORT

! parallel: enable notify instead of global barrier
#undef FMM_NOTIFY

! FMM internal preprocessors (do not change)
#undef FMM_PASS3IJKB
#undef FMM_IBOXUPD3
#undef FMM_GETDIST
#undef FMM_UNIFORMGRID
#undef FMM_LINEAR_POTENTIAL_CUBIC_EXTENSION
#undef FMM_Z_EDITING
#undef FMM_NO_INTERPOLATION_IN_NON_GRIDBOX
#undef FMM_DAMPING_SHQEN
